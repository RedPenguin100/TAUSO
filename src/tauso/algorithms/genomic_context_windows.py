import pandas as pd
from typing import Dict, Iterable
from notebooks.consts import SEQUENCE, CANONICAL_GENE
from tauso.genome.TranscriptMapper import GeneCoordinateMapper
from tauso.util import _to_str_seq

def add_external_mrna_and_context_columns(
        df: pd.DataFrame,
        mapper: GeneCoordinateMapper,
        gene_registry: Dict[str, Dict[str, str]],
        *,
        flank_sizes_premrna: Iterable[int] = (20, 30, 40, 50, 60, 70),
        flank_sizes_cds: Iterable[int] = (20, 30, 40, 50, 60, 70),
) -> pd.DataFrame:
    df = df.copy()

    SENSE_START, IN_CODING = "sense_start", "in_coding_region"
    df[IN_CODING] = False

    flank_sizes_premrna = list(flank_sizes_premrna)
    flank_sizes_cds = list(flank_sizes_cds)

    for fs in flank_sizes_premrna: df[f"flank_sequence_{fs}"] = ""
    for fs in flank_sizes_cds: df[f"local_coding_region_around_ASO_{fs}"] = ""

    for ridx, row in df.iterrows():
        idx = row.get(SENSE_START)
        if pd.isna(idx) or idx == -1: continue
        idx = int(idx)

        gene = row[CANONICAL_GENE]
        seqs = gene_registry.get(gene)
        if not seqs: continue

        pre_mrna = seqs['pre_mrna_sequence']
        cds_seq = seqs['cds_sequence']  # Now guaranteed to be pure CDS (ATG start)
        aso_len = len(_to_str_seq(row[SEQUENCE]))

        # 2. Pre-mRNA Context (Unchanged)
        for fs in flank_sizes_premrna:
            start, end = max(0, idx - fs), min(len(pre_mrna), idx + aso_len + fs)
            df.at[ridx, f"flank_sequence_{fs}"] = pre_mrna[start:end]

        # 3. Map to Spliced Index & CDS Context
        # spliced_idx is now relative to ATG (index 0)
        spliced_idx = mapper.get_spliced_index(gene, idx)

        if spliced_idx is not None and 0 <= spliced_idx < len(cds_seq):
            df.at[ridx, IN_CODING] = True

            for fs in flank_sizes_cds:
                # Raw window boundaries
                raw_start = max(0, spliced_idx - fs)
                raw_end = min(len(cds_seq), spliced_idx + aso_len + fs)

                # --- FIX: FRAME SNAPPING ---
                # Ensure the window starts at a codon boundary (multiple of 3)
                # so that calc_CAI reads correct triplets.
                frame_correction = raw_start % 3
                aligned_start = raw_start - frame_correction

                # We extend the window slightly to the left to maintain alignment
                df.at[ridx, f"local_coding_region_around_ASO_{fs}"] = cds_seq[aligned_start:raw_end]

    # 4. Final Flags
    for fs in flank_sizes_cds:
        df[f"region_is_local_{fs}"] = (df[f"local_coding_region_around_ASO_{fs}"].str.len() > 0).astype(int)

    return df