from typing import Dict, Iterable

import pandas as pd

from ..data.consts import *
from ..data.consts import SENSE_START
from ..genome.TranscriptMapper import GeneCoordinateMapper

# Assuming SENSE_START, CANONICAL_GENE, SEQUENCE, etc. are imported


def add_external_mrna_and_context_columns(
    df: pd.DataFrame,
    mapper: "GeneCoordinateMapper",  # Quoted to avoid NameError if not imported here
    gene_registry: Dict[str, Dict[str, str]],
    *,
    flank_sizes_premrna: Iterable[int] = (20, 30, 40, 50, 60, 70),
    flank_sizes_cds: Iterable[int] = (20, 30, 40, 50, 60, 70),
) -> pd.DataFrame:
    df = df.copy()

    # --- LOUD FAILURE 1: Check Inputs ---
    required_cols = [SENSE_START, CANONICAL_GENE, SEQUENCE]
    missing = [c for c in required_cols if c not in df.columns and c != SEQUENCE]
    if missing:
        raise ValueError(
            f"❌ CRITICAL ERROR: Input DataFrame is missing required columns: {missing}"
        )

    if not gene_registry:
        raise ValueError(
            "❌ CRITICAL ERROR: 'gene_registry' is empty! Cannot lookup sequences."
        )

    IN_CODING = "in_coding_region"
    flank_sizes_premrna = list(flank_sizes_premrna)
    flank_sizes_cds = list(flank_sizes_cds)

    # --- PERFORMANCE FIX: Pre-allocate standard Python lists ---
    n_rows = len(df)
    in_coding_list = [False] * n_rows

    # Dictionaries holding lists of strings for each flank size
    premrna_cols = {fs: [""] * n_rows for fs in flank_sizes_premrna}
    cds_cols = {fs: [""] * n_rows for fs in flank_sizes_cds}

    populated_count = 0

    # Extracting columns to numpy/lists and using zip is vastly faster than iterrows()
    sense_starts = df[SENSE_START].values
    genes = df[CANONICAL_GENE].values
    sequences = df[SEQUENCE].values

    for i, (idx, gene, aso_seq) in enumerate(zip(sense_starts, genes, sequences)):
        # 1. Validate Start Index
        if pd.isna(idx) or idx == -1:
            continue

        idx = int(idx)

        # 2. Validate Gene
        seqs = gene_registry.get(gene)
        if not seqs:
            continue

        pre_mrna = seqs.get("pre_mrna_sequence", "")
        cds_seq = seqs.get("cds_sequence", "")

        if not pre_mrna:
            continue

        aso_len = len(str(aso_seq))  # Safe access

        # --- POPULATE ---
        populated_count += 1

        # Pre-mRNA Context
        for fs in flank_sizes_premrna:
            start, end = max(0, idx - fs), min(len(pre_mrna), idx + aso_len + fs)
            premrna_cols[fs][i] = pre_mrna[start:end]

        # CDS Context
        spliced_idx = mapper.get_spliced_index(gene, idx)

        if spliced_idx is not None and 0 <= spliced_idx < len(cds_seq):
            in_coding_list[i] = True
            for fs in flank_sizes_cds:
                raw_start = max(0, spliced_idx - fs)
                raw_end = min(len(cds_seq), spliced_idx + aso_len + fs)

                # Frame alignment
                frame_correction = raw_start % 3
                aligned_start = raw_start - frame_correction

                cds_cols[fs][i] = cds_seq[aligned_start:raw_end]

    # --- LOUD FAILURE 3: Check Result ---
    if populated_count == 0:
        raise RuntimeError(
            "❌ CRITICAL FAILURE: Function finished but 0 sequences were populated.\n"
            "   Likely Causes:\n"
            "   1. All 'sense_start' values are NaN or -1.\n"
            "   2. No gene names in DF match keys in 'gene_registry'.\n"
            "   3. 'gene_registry' values are missing 'pre_mrna_sequence'."
        )

    # --- PERFORMANCE FIX: Assign lists back to DataFrame in bulk ---
    df[IN_CODING] = in_coding_list
    for fs in flank_sizes_premrna:
        df[f"flank_sequence_{fs}"] = premrna_cols[fs]
    for fs in flank_sizes_cds:
        df[f"local_coding_region_around_ASO_{fs}"] = cds_cols[fs]

    # 4. Final Flags
    for fs in flank_sizes_cds:
        df[f"region_is_local_{fs}"] = (
            df[f"local_coding_region_around_ASO_{fs}"].str.len() > 0
        ).astype(int)

    return df
