import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

from ..data.consts import CANONICAL_GENE, CELL_LINE, standardize_cell_line_name
from ..data.data import get_data_dir
from ..features.codon_usage.cai import calc_CAI
from ..features.codon_usage.enc import compute_ENC
from ..features.codon_usage.tai import calc_tAI, tai_weights

logger = logging.getLogger(__name__)


def populate_tai(df: pd.DataFrame, cds_windows: list, registry: dict) -> tuple[pd.DataFrame, list[str]]:
    """
    Calculates local and global tAI scores.
    Raises KeyError if required local context columns are missing.
    """
    weights = tai_weights("hm")
    feature_names = []

    # 1. Local tAI Scores
    logger.info("Calculating local tAI scores for CDS windows...")
    for flank in cds_windows:
        local_col = f"local_coding_region_around_ASO_{flank}"

        # Strict check: Raise error if column is missing
        if local_col not in df.columns:
            raise KeyError(
                f"Missing required column for tAI calculation: '{local_col}'. Please run add_external_mrna_and_context_columns"
            )

        feat_name = f"tAI_score_{flank}_CDS"
        df[feat_name] = df[local_col].apply(lambda seq: calc_tAI(seq, weights))
        feature_names.append(feat_name)

    # 2. Global tAI Scores (Registry Lookup)
    logger.info("Pre-calculating Global tAI for %d unique genes...", len(registry))
    gene_tai_lookup = {
        gene: calc_tAI(data.get("cds_sequence"), weights) if data.get("cds_sequence") else np.nan
        for gene, data in registry.items()
    }

    logger.info("Mapping Global tAI scores to dataframe...")
    df["tAI_score_global_CDS"] = df[CANONICAL_GENE].map(gene_tai_lookup)
    feature_names.append("tAI_score_global_CDS")

    return df, feature_names


def populate_enc(
    df: pd.DataFrame,
    cds_windows: list,
    registry: dict,
    n_jobs: int = 1,
) -> tuple[pd.DataFrame, list[str]]:
    """
    Calculates local and global ENC (Effective Number of Codons) scores.

    Args:
        df: Input dataframe.
        cds_windows: List of window sizes.
        registry: Dictionary of gene info.
        n_jobs: Number of CPU nodes to use. Defaults to 1 (standard .apply).

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
            pandarallel.initialize(nb_workers=n_jobs)
        except ImportError:
            logger.warning("'pandarallel' not found. Falling back to single core.")
            use_parallel = False

    # 2. Local ENC Scores
    logger.info("Calculating local ENC scores for CDS windows...")

    for flank in cds_windows:
        local_col = f"local_coding_region_around_ASO_{flank}"
        enc_col = f"ENC_score_{flank}_CDS"

        if local_col not in df.columns:
            raise KeyError(
                f"Missing required column for ENC calculation: '{local_col}'.  Please run add_external_mrna_and_context_columns"
            )

        if use_parallel:
            df[enc_col] = df[local_col].parallel_apply(compute_ENC)
        else:
            df[enc_col] = df[local_col].apply(compute_ENC)

        feature_names.append(enc_col)

    # 3. Global ENC Scores (Registry Lookup)
    logger.info("Pre-calculating Global ENC for %d unique genes...", len(registry))

    gene_enc_lookup = {
        gene: compute_ENC(data.get("cds_sequence")) if data.get("cds_sequence") else np.nan
        for gene, data in registry.items()
    }

    logger.info("Mapping Global ENC scores to dataframe...")

    df["ENC_score_global_CDS"] = df[CANONICAL_GENE].map(gene_enc_lookup)
    feature_names.append("ENC_score_global_CDS")

    logger.info("Global ENC Calculation Complete.")

    return df, feature_names


def populate_cai(
    df: pd.DataFrame,
    cds_windows: list,
    registry: dict,
    n_jobs: int = 1,
) -> tuple[pd.DataFrame, list[str]]:
    """
    Calculates context-specific CAI scores based on cell line weights.
    Loads 'cai_weights.json' internally.

    Args:
        df: Input dataframe.
        cds_windows: List of window sizes.
        registry: Dictionary of gene info (ref_registry).
        n_jobs: Number of CPU nodes (defaults to 1).

    Returns:
        tuple: (Updated DataFrame, List of new feature names)
    """
    feature_names = []

    # 1. Load Centralized Weights
    weights_path = Path(get_data_dir()) / "cai_weights.json"
    if not weights_path.exists():
        raise FileNotFoundError(f"CAI weights file not found at: {weights_path}")

    with open(weights_path, "r") as f:
        weight_map = json.load(f)

    logger.info("Loaded CAI profiles for: %s", list(weight_map.keys()))

    # 2. Setup Parallelization
    use_parallel = n_jobs > 1
    if use_parallel:
        try:
            from pandarallel import pandarallel

            pandarallel.initialize(nb_workers=n_jobs)
        except ImportError:
            logger.warning("'pandarallel' not found. Falling back to single core.")
            use_parallel = False

    # 3. Warn once per unknown cell line before any parallel work.
    # pandarallel pickles the closure into worker processes, so per-row
    # deduplication sets would fire once per worker instead of once total.
    _DEPMAP_NO_EXPRESSION = {"SW872", "HK2"}
    unknown_cell_lines = (
        set(df[CELL_LINE].map(standardize_cell_line_name).dropna()) - set(weight_map.keys())
    )
    for _cl in sorted(unknown_cell_lines):
        if _cl in _DEPMAP_NO_EXPRESSION:
            logger.warning(
                "Cell line '%s' not found in weights (found in DepMap but lacks expression data). "
                "Using 'Generic' profile.",
                _cl,
            )
        else:
            logger.warning("Cell line '%s' not found in weights. Using 'Generic' profile.", _cl)

    def _get_row_cai(row, sequence):
        """Internal helper to calculate CAI for a single row."""
        if pd.isna(sequence) or not str(sequence).strip():
            return np.nan

        cell_line_raw = row.get(CELL_LINE)
        cell_line_name = standardize_cell_line_name(cell_line_raw)

        weights = weight_map.get(cell_line_name) or weight_map.get("Generic")

        if not weights:
            return np.nan

        return calc_CAI(str(sequence), weights)

    # 4. Local CAI Windows
    logger.info("Applying context-specific CAI scores to windows...")

    for flank in cds_windows:
        local_col = f"local_coding_region_around_ASO_{flank}"
        cai_col = f"CAI_score_{flank}_CDS"

        if local_col not in df.columns:
            raise KeyError(
                f"Missing required column for cAI calculation: '{local_col}'."
                f" Please run add_external_mrna_and_context_columns"
            )

        # Define the lambda for this specific column, binding local_col immediately
        # Note: axis=1 is required because we need both sequence (col) and cell_line (row)
        op = lambda row, local_col=local_col: _get_row_cai(row, row[local_col])

        if use_parallel:
            df[cai_col] = df.parallel_apply(op, axis=1)
        else:
            df[cai_col] = df.apply(op, axis=1)

        feature_names.append(cai_col)

    # 5. Global CAI (Registry Lookup)
    logger.info("Calculating Global CDS CAI Scores...")

    def _get_global_cai(row):
        gene = row.get(CANONICAL_GENE)
        entry = registry.get(gene)

        if not entry or not entry.get("cds_sequence"):
            return np.nan

        return _get_row_cai(row, entry["cds_sequence"])

    if use_parallel:
        df["CAI_score_global_CDS"] = df.parallel_apply(_get_global_cai, axis=1)
    else:
        df["CAI_score_global_CDS"] = df.apply(_get_global_cai, axis=1)

    feature_names.append("CAI_score_global_CDS")

    logger.info("CAI Calculation Complete.")

    return df, feature_names
