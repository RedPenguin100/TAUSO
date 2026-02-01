import numpy as np
import pandas as pd

from ..data.consts import CANONICAL_GENE
from ..features.codon_usage.tai import calc_tAI, tai_weights
from ..features.codon_usage.enc import compute_ENC


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
            raise KeyError(f"Missing required column for ENC calculation: '{local_col}'")

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
