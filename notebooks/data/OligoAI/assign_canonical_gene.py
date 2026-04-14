import pandas as pd
import warnings

from notebooks.consts import ORIGINAL_OLIGO_CSV, ORIGINAL_OLIGO_CSV_WITH_CANONICAL
from tauso.data.data import load_db, get_paths
from tauso.off_target.search import find_all_gene_off_targets, get_bowtie_index_base
from tauso.data.consts import CANONICAL_GENE

import pandas as pd
import warnings


import pandas as pd
import warnings
from collections import Counter
import concurrent.futures
from tqdm import tqdm  # Highly recommended for tracking parallel progress



import pandas as pd
import warnings
import uuid
import os
from collections import Counter

import subprocess
import pandas as pd
from collections import defaultdict

import subprocess
import os
import pandas as pd
from collections import defaultdict
import logging

logger = logging.getLogger(__name__)


def find_all_gene_off_targets_BULK(fasta_path, genome="GRCh38", threads=16, max_mismatches=0):
    """
    Main bulk entry point. Replaces your old find_all_gene_off_targets.
    """
    # 1. Run the bulk Bowtie search
    hits_list = run_bowtie_search_bulk(fasta_path, genome=genome, max_mismatches=max_mismatches, threads=threads)

    # 2. Annotate the hits in bulk and get the mapping dictionary
    seq_to_genes = annotate_hits_bulk(hits_list, genome=genome)

    return seq_to_genes


def run_bowtie_search_bulk(fasta_path, genome="GRCh38", max_mismatches=0, threads=16):
    """
    Runs Bowtie 1 alignment on a whole FASTA file.
    """
    index_base = get_bowtie_index_base(genome=genome)
    sam_output = fasta_path.replace(".fasta", ".sam")

    cmd = [
        "bowtie",
        "-v",
        str(max_mismatches),
        "-a",  # Report all valid alignments
        "-S",
        "--sam-nohead",
        "-p",
        str(threads),  # MULTITHREADING ENABLED
        "-f",  # FASTA INPUT (changed from -c)
        "-x",
        index_base,
        fasta_path,  # Input file
        sam_output,  # Output file
    ]

    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Bowtie bulk search failed: {e.stderr}")
        return []

    hits = []

    # Parse the massive SAM output line by line
    with open(sam_output, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("@"):
                continue

            parts = line.split("\t")
            flag = int(parts[1])
            if flag & 4:
                continue  # Unmapped

            # Because we wrote the FASTA with >SEQUENCE, parts[0] is the sequence itself
            sequence_id = parts[0]
            chrom = parts[2]
            start_pos = int(parts[3]) - 1

            hits.append(
                {"sequence": sequence_id, "chrom": chrom, "start": start_pos, "end": start_pos + len(sequence_id)}
            )

    if os.path.exists(sam_output):
        os.remove(sam_output)

    return hits


import pyranges as pr
import pandas as pd
import time

def load_pyranges_db(gtf_path):
    print("Loading GTF into RAM... (This takes ~45s but only happens once)")
    gr = pr.read_gtf(gtf_path)

    # Tell pyranges to mathematically infer introns from the exons
    if "intron" not in gr.Feature.unique():
        print("Inferring missing introns from transcript coordinates...")
        gr_introns = gr.features.introns()
        gr = pr.concat([gr, gr_introns])

    return gr

import pyranges as pr
import pandas as pd
import time

def annotate_hits_bulk(hits_list, genome):
    """
    Lightning-fast pyranges intersection replacing the SQLite loop.
    Pure spatial intersection: returns ALL genes an ASO touches,
    ignoring feature types, introns, or biotypes.
    """
    if not hits_list:
        return {}

    # Assumes get_paths is defined globally in your script
    gr_genome = load_pyranges_db(get_paths(genome)['gtf'])

    start_time = time.time()

    # 1. Convert your raw Bowtie hits to a Pandas DataFrame
    df_hits = pd.DataFrame(hits_list)

    # Pyranges explicitly requires these exact column names (capitalized)
    df_hits = df_hits.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})

    # Convert to a PyRanges object
    gr_hits = pr.PyRanges(df_hits)

    print(f"Intersecting {len(df_hits)} hits against the genome...")

    # 2. THE C++ INTERSECTION (Takes ~2 seconds)
    # Added apply_strand_suffix=False to silence the strand warning
    intersected = gr_hits.join(gr_genome, apply_strand_suffix=False)
    df_res = intersected.df

    if df_res.empty:
        return {}

    # 3. Resolve the gene name based on the fallback chain
    if "gene_name" in df_res.columns:
        df_res["resolved_gene_name"] = df_res["gene_name"]
    elif "Name" in df_res.columns:
        df_res["resolved_gene_name"] = df_res["Name"]
    elif "gene_id" in df_res.columns:
        df_res["resolved_gene_name"] = df_res["gene_id"]
    else:
        df_res["resolved_gene_name"] = None

    # Drop anything that didn't map to a name
    df_res = df_res.dropna(subset=["resolved_gene_name"])

    # 4. AGGREGATE BACK TO DICTIONARY (No sorting or dropping duplicates needed)
    # Group by the ASO sequence and collect all unique gene names it hit
    seq_to_genes_series = df_res.groupby("sequence")["resolved_gene_name"].unique()

    # Convert pandas Series of arrays to a standard Python dictionary of lists
    seq_to_genes = {seq: list(genes) for seq, genes in seq_to_genes_series.items()}

    elapsed = time.time() - start_time
    print(f"Annotation complete in {elapsed:.2f} seconds!")

    return seq_to_genes
#
# def annotate_hits_bulk(hits_list, genome="GRCh38"):
#     """
#     Annotates hits with a high-speed cache and returns the dictionary mapping needed.
#     """
#     if not hits_list:
#         return {}
#
#     db = load_db(genome=genome)
#     seq_to_genes = defaultdict(set)
#
#     # CACHE: Drastically reduces SQLite db.region() overhead
#     region_cache = {}
#
#     print(f"Annotating {len(hits_list)} total genome hits...")
#
#     # Wrap the iteration in tqdm for a progress bar
#     for hit in tqdm(hits_list, desc="Annotating hits"):
#         seq = hit["sequence"]
#         chrom = hit["chrom"]
#         start = hit["start"]
#         end = hit["end"]
#
#         region_key = (chrom, start, end)
#
#         # If we've already looked up this exact coordinate, use the cache
#         if region_key in region_cache:
#             gene_name = region_cache[region_key]
#         else:
#             try:
#                 features = list(db.region(region=(chrom, start, end)))
#             except Exception:
#                 features = []
#
#             gene_name = None
#             current_priority = 0
#
#             # Priority logic
#             for feat in features:
#                 f_type = feat.featuretype
#
#                 if f_type == "exon":
#                     priority = 4
#                 elif f_type == "CDS":
#                     priority = 4
#                 elif f_type == "intron":
#                     priority = 2
#                 elif f_type == "gene":
#                     priority = 1
#                 else:
#                     priority = 0
#
#                 if priority > current_priority:
#                     current_priority = priority
#                     gene_name = feat.attributes.get(
#                         "gene_name",
#                         feat.attributes.get("Name", feat.attributes.get("gene_id", [None])),
#                     )[0]
#
#             # Save to cache for the next time this coordinate appears
#             region_cache[region_key] = gene_name
#
#         # Map the found gene to the ASO sequence
#         if gene_name:
#             seq_to_genes[seq].add(gene_name)
#
#     # Convert the sets to lists before returning to match the main script's expectations
#     return {seq: list(genes) for seq, genes in seq_to_genes.items()}


import pandas as pd
import warnings
import uuid
import os
from collections import Counter


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

# 1. Worker function MUST be at the top level (global scope) for multiprocessing
def _process_gene_group(args):
    target_gene, seqs, genome = args

    all_found_genes = set()
    hit_counts = Counter()

    for aso_seq in seqs:
        dna = aso_seq.upper().replace("U", "T")
        # Assuming find_all_gene_off_targets is globally available
        hits = find_all_gene_off_targets(dna, genome=genome, max_mismatches=0)

        if not hits.empty:
            genes = hits["gene_name"].dropna().unique()
            all_found_genes.update(genes)

            # Tally up every gene found by this specific ASO
            for g in genes:
                hit_counts[g] += 1

    # Return the target, the aggregated set, the raw counts, and total sequences tested
    return target_gene, list(all_found_genes), dict(hit_counts), len(seqs)


def process_and_assign_genes(
    input_csv=ORIGINAL_OLIGO_CSV,
    output_csv=ORIGINAL_OLIGO_CSV_WITH_CANONICAL,
    genome="GRCh38",
    sequence_col="rna_sequence" # <-- ADDED: specify your ASO column name here
):
    df = pd.read_csv(input_csv)

    # Changed 'rna_context' to your sequence column
    missing_cols = {sequence_col, "target_gene"} - set(df.columns)
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    # Group by target_gene and grab the first ASO sequence instead
    gene_to_seq = df.dropna(subset=["target_gene", sequence_col]).groupby("target_gene")[sequence_col].first()

    print(f"Resolving {len(gene_to_seq)} unique target genes")

    gene_to_canonical = {}
    ambiguous_genes = {}

    for i, (target_gene, aso_seq) in enumerate(gene_to_seq.items(), start=1):
        # aso_seq might already be DNA, but replacing U->T is a safe fallback
        dna = aso_seq.upper().replace("U", "T")
        hits = find_all_gene_off_targets(dna, genome=genome, max_mismatches=0)
        genes = hits["gene_name"].dropna().unique().tolist() if not hits.empty else []

        if len(genes) == 0:
            gene_to_canonical[target_gene] = pd.NA
        elif len(genes) == 1:
            gene_to_canonical[target_gene] = genes[0]
        else:
            gene_to_canonical[target_gene] = ";".join(genes)
            ambiguous_genes[target_gene] = genes

        if i % 50 == 0 or i == len(gene_to_seq):
            print(f"[{i}/{len(gene_to_seq)}] genes processed")

    df[CANONICAL_GENE] = df["target_gene"].map(gene_to_canonical)
    df.to_csv(output_csv, index=False)
    print(f"Saved annotated dataframe → {output_csv}")

    if ambiguous_genes:
        warnings.warn(f"{len(ambiguous_genes)} target genes mapped to multiple canonical genes")
    else:
        print("No ambiguous gene mappings found 🎉")

    return df, ambiguous_genes


if __name__ == "__main__":
    # This block executes only when run from the command line (e.g., in CI)
    process_and_assign_genes()