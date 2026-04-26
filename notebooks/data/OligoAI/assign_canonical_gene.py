import ast
import logging
import os
import uuid
from collections import Counter

import pandas as pd

from notebooks.consts import ORIGINAL_OLIGO_CSV, ORIGINAL_OLIGO_CSV_WITH_CANONICAL
from tauso.data.consts import CANONICAL_GENE
from tauso.off_target.search import find_all_gene_off_targets_BULK

logger = logging.getLogger(__name__)


def process_and_assign_genes_bulk(
    input_csv=ORIGINAL_OLIGO_CSV,
    output_csv=ORIGINAL_OLIGO_CSV_WITH_CANONICAL,
    genome="GRCh38",
    sequence_col="rna_sequence",
    threads=16,
):
    df = pd.read_csv(input_csv)

    missing_cols = {sequence_col, "target_gene"} - set(df.columns)
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    # 1. Extract and write bulk FASTA
    all_unique_seqs = df[sequence_col].dropna().unique()
    print(f"Extracted {len(all_unique_seqs)} unique ASO sequences.")

    run_id = uuid.uuid4().hex
    temp_fasta = f"bulk_input_{run_id}.fasta"

    with open(temp_fasta, "w") as f:
        for seq in all_unique_seqs:
            dna = seq.upper().replace("U", "T")
            f.write(f">{dna}\n{dna}\n")

    # 2. Run bulk alignment
    print("Running bulk alignment and annotation...")
    sequence_to_genes_map = find_all_gene_off_targets_BULK(temp_fasta, genome, threads)

    if os.path.exists(temp_fasta):
        os.remove(temp_fasta)

    # 3. Map results and calculate stats
    print("Mapping results back to dataframe and calculating stats...")
    gene_to_seqs = df.dropna(subset=["target_gene", sequence_col]).groupby("target_gene")[sequence_col].unique()

    gene_to_canonical = {}
    ambiguous_genes = {}
    gene_stats = {}

    for target_gene, seqs in gene_to_seqs.items():
        all_found_genes = set()
        hit_counts = Counter()

        for aso_seq in seqs:
            dna = aso_seq.upper().replace("U", "T")
            genes = sequence_to_genes_map.get(dna, [])

            if genes:
                all_found_genes.update(genes)
                for g in genes:
                    hit_counts[g] += 1

        # --- Populate Statistics ---
        total_seqs = len(seqs)
        if hit_counts:
            most_common_gene = max(hit_counts, key=hit_counts.get)
            most_common_count = hit_counts[most_common_gene]
        else:
            most_common_gene = pd.NA
            most_common_count = 0

        gene_stats[target_gene] = {
            "total_asos_tested": total_seqs,
            "most_popular_canonical": most_common_gene,
            "popular_hit_count": most_common_count,
            "match_rate": f"{most_common_count}/{total_seqs}",
            "all_hit_frequencies": dict(hit_counts),
        }

    df.to_csv(output_csv, index=False)
    print(f"Saved annotated dataframe → {output_csv}")

    # 4. Convert stats to DataFrame and export
    df_stats = pd.DataFrame.from_dict(gene_stats, orient="index").reset_index()
    df_stats = df_stats.rename(columns={"index": "original_target_gene"})

    stats_csv = str(output_csv).replace(".csv", "_STATS.csv")
    df_stats.to_csv(stats_csv, index=False)
    print(f"Saved detailed stats dataframe → {stats_csv}")

    return df, ambiguous_genes, df_stats


# The complete mapping dictionary for the 24 mismatched groups
MANUAL_CANONICAL_MAPPING = {
    # Original is wrong (compilation error)
    "SFRS10": "TRA2B",  # Outdated alias
    "SENP2": "SUMO2",  # Patent targeted SUMO2 but compiler confused it with enzyme that cuts it, SENP2
    "SRC": "CSK",  # Patent targeted C-SRC tyrosine kinase which is CSK but GPT stopped at SRC
    "NRF1": "NKRF",  # Patent targeted "NRF" (NF-kappa-B-repressing factor, now NKRF) but GPT confused with NRF1
    "BRCA2": "N4BP2L2",  # Patent targeted BRCA2 region transcription unit CG005 which later renamed to N4BP2L2, GPT ignored the later fact.
    "MAPK6": "MAPK12",  # Patent targeted ERK6, historic name for MAPK12, but compiler confused with ERK3 which was mentioned in the patent
    "NFKBIL1": "TONSL",  # Patent targeted IKB-R whihc is the old name for NFKBIL2 and officially TONSL, GPT guessed wrong
    "PTPRJ": "PTPRF",  # Patent target LAR which is alias for PTPRF, GPT mixed PTP receptor with PTPRJ
    "THRAP3": "ZNHIT3",  # Patent targeted TRIP3, which is an alias for ZNHIT3, but mistakenly transformed TRIP3 to THRAP3.
    # Original is wrong (patent / mislabel)
    "UBE3A": "SNHG14",  # Patent *discusses* UBE3A regulation, but targets SNHG14 completely different gene UNKL
    "RASD1": "RASD2",  # Patent probably annotated incorrectly
    # The most popular disagrees with the original, we choose the original
    "ANGPT2": "ANGPT2",
    "ANGPTL3": "ANGPTL3",
    "BAG5": "BAG5",
    "FNTB": "FNTB",
    "IGF2": "IGF2",
    "IL18": "IL18",
    "MMP1": "MMP1",
    "MRPS16": "MRPS16",
    "PLN": "PLN",
    "S1PR2": "S1PR2",
    "STAT5B": "STAT5B",
    # Patent errors, unrecoverable
    "RHO": pd.NA,  # Patent targeting mutation, we don't mutate the genome.
    "NAIP": pd.NA,  # Patent probably contains an error (Old, genome wasn't annotated properly). ASOs target
}


def refine_and_clean_canonical_genes(df, gene_stats, output_csv, canonical_mapping=MANUAL_CANONICAL_MAPPING):
    """
    Analyzes hit statistics, prints discrepancy reports, and assigns
    the final CANONICAL_GENE using original targets and manual overrides.
    """
    # If the most popular candidate agrees with the original candidate we can safely continue.
    filtered_df = gene_stats[gene_stats["original_target_gene"] != gene_stats["most_popular_canonical"]].copy()

    # 1. Safely parse the dictionary column
    def parse_dict(val):
        if isinstance(val, str):
            try:
                return ast.literal_eval(val)
            except:
                return {}
        return val if isinstance(val, dict) else {}

    filtered_df["all_hit_frequencies"] = filtered_df["all_hit_frequencies"].apply(parse_dict)

    # 2. Extract the hit count for the original target gene
    def get_original_count(row):
        original = row["original_target_gene"]
        freq_dict = row["all_hit_frequencies"]
        return freq_dict.get(original, 0)

    filtered_df["original_hit_count"] = filtered_df.apply(get_original_count, axis=1)

    # 3. Use a 50% threshold to ensure we don't keep noise (like NAIP's 3 hits)
    threshold = filtered_df["total_asos_tested"] * 0.5
    present_originals = filtered_df[filtered_df["original_hit_count"] >= threshold]
    missing_originals = filtered_df[filtered_df["original_hit_count"] < threshold]

    # 4. Print the MISSING ORIGINALS
    print(f"--- MISSING ORIGINALS ({len(missing_originals)} groups) ---")
    print("These did not map to the original name at all (or missed the 50% threshold):\n")
    missing_cols = [
        "original_target_gene",
        "original_hit_count",
        "total_asos_tested",
        "most_popular_canonical",
        "popular_hit_count",
        "match_rate",
    ]
    if not missing_originals.empty:
        print(missing_originals[missing_cols].to_string(index=False))

    # 5. Print the RECOVERABLE ORIGINALS
    print(f"\n--- RECOVERABLE ORIGINALS ({len(present_originals)} groups) ---")
    print("These exist in the frequencies and can be reverted to the original name:\n")
    recoverable_cols = [
        "original_target_gene",
        "original_hit_count",
        "total_asos_tested",
        "most_popular_canonical",
        "popular_hit_count",
    ]
    if not present_originals.empty:
        print(present_originals[recoverable_cols].to_string(index=False))

    # --- APPLYING CORRECTIONS ---
    print("\n--- APPLYING CORRECTIONS & CLEANING DATASET ---")
    initial_count = len(df)

    # 1. First, drop ASOs that don't even have a target_gene to begin with
    df_no_empty_genes = df.dropna(subset=["target_gene"]).copy()
    dropped_empty = initial_count - len(df_no_empty_genes)
    print(f"Dropped {dropped_empty} ASOs missing an initial 'target_gene'.")

    # 2. Apply the map using .replace() to preserve our pd.NA traps
    df_no_empty_genes[CANONICAL_GENE] = df_no_empty_genes["target_gene"].replace(canonical_mapping)

    # 3. Purge the unrecoverable data (the ones we mapped to pd.NA like RHO, NAIP)
    df_clean = df_no_empty_genes.dropna(subset=[CANONICAL_GENE]).copy()
    dropped_unrecoverable = len(df_no_empty_genes) - len(df_clean)

    print(f"Dataset cleaned! Purged {dropped_unrecoverable} ASOs due to unrecoverable errors (e.g., RHO, NAIP).")
    print(f"Final dataset contains {len(df_clean)} ASOs ready for TAUSO.")

    # 4. Save to file
    df_clean.to_csv(output_csv, index=False)
    print(f"Successfully saved to {output_csv}")

    return df_clean


if __name__ == "__main__":
    # 1. Run the heavy alignment computation to generate stats
    df, ambiguous_genes, gene_stats = process_and_assign_genes_bulk(threads=16)

    # 2. Use those stats to analyze discrepancies and assign the final canonical genes
    df_clean = refine_and_clean_canonical_genes(
        df=df, gene_stats=gene_stats, output_csv=ORIGINAL_OLIGO_CSV_WITH_CANONICAL
    )

    # Optional: Quick check for DMPK if you still want to inspect it directly
