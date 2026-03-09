import json
from pathlib import Path

import numpy as np
import pandas as pd

from ..data.consts import CANONICAL_GENE, CELL_LINE, standardize_cell_line_name
from ..data.data import get_data_dir
from ..features.codon_usage.cai import calc_CAI
from ..features.codon_usage.enc import compute_ENC
from ..features.codon_usage.tai import calc_tAI, tai_weights


def populate_tai(df: pd.DataFrame, cds_windows: list, registry: dict) -> tuple[pd.DataFrame, list[str]]:
    """
    Calculates local and global tAI scores.
    Raises KeyError if required local context columns are missing.
    """
    weights = tai_weights("hm")
    feature_names = []

    # 1. Local tAI Scores
    print("Calculating local tAI scores for CDS windows...")
    for flank in cds_windows:
        local_col = f"local_coding_region_around_ASO_{flank}"

        # Strict check: Raise error if column is missing
        if local_col not in df.columns:
            raise KeyError(
                f"Missing required column for tAI calculation: '{local_col}'. Please run add_external_mrna_and_context_columns")

        feat_name = f"tAI_score_{flank}_CDS"
        df[feat_name] = df[local_col].apply(lambda seq: calc_tAI(seq, weights))
        feature_names.append(feat_name)

    # 2. Global tAI Scores (Registry Lookup)
    print(f"Pre-calculating Global tAI for {len(registry)} unique genes...")
    gene_tai_lookup = {
        gene: calc_tAI(data.get('cds_sequence'), weights) if data.get('cds_sequence') else np.nan
        for gene, data in registry.items()
    }

    print("Mapping Global tAI scores to dataframe...")
    df["tAI_score_global_CDS"] = df[CANONICAL_GENE].map(gene_tai_lookup)
    feature_names.append("tAI_score_global_CDS")

    return df, feature_names


def populate_enc(df: pd.DataFrame, cds_windows: list, registry: dict,
                 n_jobs: int = 1, verbose: bool = False) -> tuple[pd.DataFrame, list[str]]:
    """
    Calculates local and global ENC (Effective Number of Codons) scores.

    Args:
        df: Input dataframe.
        cds_windows: List of window sizes.
        registry: Dictionary of gene info.
        n_jobs: Number of CPU nodes to use. Defaults to 1 (standard .apply).
        verbose: If True, prints progress messages.

    Returns:
        tuple: (Updated DataFrame, List of new feature names)
    """
    feature_names = []

    # 1. Setup Parallelization
    use_parallel = n_jobs > 1
    if use_parallel:
        try:
            from pandarallel import pandarallel
            # Initialize strictly if using parallelism
            pandarallel.initialize(nb_workers=n_jobs, progress_bar=verbose, verbose=0)
        except ImportError:
            if verbose: print("Warning: 'pandarallel' not found. Falling back to single core.")
            use_parallel = False

    # 2. Local ENC Scores
    if verbose: print("Calculating local ENC scores for CDS windows...")

    for flank in cds_windows:
        local_col = f"local_coding_region_around_ASO_{flank}"
        enc_col = f"ENC_score_{flank}_CDS"

        if local_col not in df.columns:
            raise KeyError(
                f"Missing required column for ENC calculation: '{local_col}'.  Please run add_external_mrna_and_context_columns")

        if use_parallel:
            df[enc_col] = df[local_col].parallel_apply(compute_ENC)
        else:
            df[enc_col] = df[local_col].apply(compute_ENC)

        feature_names.append(enc_col)

    # 3. Global ENC Scores (Registry Lookup)
    if verbose: print(f"Pre-calculating Global ENC for {len(registry)} unique genes...")

    gene_enc_lookup = {
        gene: compute_ENC(data.get('cds_sequence')) if data.get('cds_sequence') else np.nan
        for gene, data in registry.items()
    }

    if verbose: print("Mapping Global ENC scores to dataframe...")

    df["ENC_score_global_CDS"] = df[CANONICAL_GENE].map(gene_enc_lookup)
    feature_names.append("ENC_score_global_CDS")

    if verbose: print("Global ENC Calculation Complete.")

    return df, feature_names


def populate_cai(df: pd.DataFrame, cds_windows: list, registry: dict,
                 n_jobs: int = 1, verbose: bool = False) -> tuple[pd.DataFrame, list[str]]:
    """
    Calculates context-specific CAI scores based on cell line weights.
    Loads 'cai_weights.json' internally.

    Args:
        df: Input dataframe.
        cds_windows: List of window sizes.
        registry: Dictionary of gene info (ref_registry).
        n_jobs: Number of CPU nodes (defaults to 1).
        verbose: Print status messages.

    Returns:
        tuple: (Updated DataFrame, List of new feature names)
    """
    feature_names = []

    # 1. Load Centralized Weights
    weights_path = Path(get_data_dir()) / "cai_weights.json"
    if not weights_path.exists():
        raise FileNotFoundError(f"CAI weights file not found at: {weights_path}")

    with open(weights_path, 'r') as f:
        weight_map = json.load(f)

    if verbose:
        print(f"Loaded CAI profiles for: {list(weight_map.keys())}")

    # 2. Setup Parallelization
    use_parallel = n_jobs > 1
    if use_parallel:
        try:
            from pandarallel import pandarallel
            pandarallel.initialize(nb_workers=n_jobs, progress_bar=verbose, verbose=0)
        except ImportError:
            if verbose: print("Warning: 'pandarallel' not found. Falling back to single core.")
            use_parallel = False

    # 3. Define Row-Wise Helper
    # We use a mutable set in the closure to track reported warnings
    reported_fallbacks = set()

    def _get_row_cai(row, sequence):
        """Internal helper to calculate CAI for a single row."""
        if pd.isna(sequence) or not str(sequence).strip():
            return np.nan

        cell_line_raw = row.get(CELL_LINE)
        cell_line_name = standardize_cell_line_name(cell_line_raw)

        # Select weights
        if cell_line_name in weight_map:
            weights = weight_map[cell_line_name]
        else:
            if verbose and cell_line_name not in reported_fallbacks:
                print(f"⚠️  Cell line '{cell_line_name}' not found in weights. Using 'Generic' profile.")
                reported_fallbacks.add(cell_line_name)
            weights = weight_map.get('Generic')

        if not weights:
            return np.nan

        return calc_CAI(str(sequence), weights)

    # 4. Local CAI Windows
    if verbose: print("Applying context-specific CAI scores to windows...")

    for flank in cds_windows:
        local_col = f"local_coding_region_around_ASO_{flank}"
        cai_col = f"CAI_score_{flank}_CDS"

        if local_col not in df.columns:
            raise KeyError(f"Missing required column for cAI calculation: '{local_col}'."
                           f" Please run add_external_mrna_and_context_columns")

        # Define the lambda for this specific column
        # Note: axis=1 is required because we need both sequence (col) and cell_line (row)
        op = lambda row: _get_row_cai(row, row[local_col])

        if use_parallel:
            df[cai_col] = df.parallel_apply(op, axis=1)
        else:
            df[cai_col] = df.apply(op, axis=1)

        feature_names.append(cai_col)

    # 5. Global CAI (Registry Lookup)
    if verbose: print("Calculating Global CDS CAI Scores...")

    def _get_global_cai(row):
        gene = row.get(CANONICAL_GENE)
        entry = registry.get(gene)

        if not entry or not entry.get('cds_sequence'):
            return np.nan

        return _get_row_cai(row, entry['cds_sequence'])

    if use_parallel:
        df["CAI_score_global_CDS"] = df.parallel_apply(_get_global_cai, axis=1)
    else:
        df["CAI_score_global_CDS"] = df.apply(_get_global_cai, axis=1)

    feature_names.append("CAI_score_global_CDS")

    if verbose: print("CAI Calculation Complete.")

    return df, feature_names
