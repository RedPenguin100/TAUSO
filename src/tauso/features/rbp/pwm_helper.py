import logging
import os

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

from ...data.consts import CELL_LINE_DEPMAP
from ...data.data import get_data_dir
from ..codon_usage.find_cai_reference import load_cell_line_gene_expression

_DATA_DIR = os.path.join(get_data_dir(), "attract")


def build_rbp_expression_matrix(
    cell_lines: list,
    pwm_db: dict,
    data_dir: str = None,
) -> pd.DataFrame:
    """
    Builds a matrix of RNA expression (TPM) for the RBPs present in the pwm_db
    across the provided list of cell lines.

    Args:
        cell_lines: List of cell line IDs (e.g., DepMap IDs).
        pwm_db: Dictionary of PWMs where keys are ATtRACT Matrix IDs.
        cell_line_col: Name to use for the cell line index in the output matrix.
        data_dir: Path to data directory. If None, uses tauso.data.data.get_data_dir().

    Returns:
        pd.DataFrame: Pivot table (Index=CellLine, Columns=Gene, Values=TPM).
    """
    # --- 1. SETUP & FAIL LOUDLY ---
    if data_dir is None:
        data_dir = get_data_dir()

    attract_csv_path = os.path.join(_DATA_DIR, "RBS_motifs_Homo_sapiens.csv")
    expression_dir = os.path.join(data_dir, "processed_expression")

    if not os.path.exists(attract_csv_path):
        raise FileNotFoundError(
            f"CRITICAL: ATtRACT metadata not found at {attract_csv_path}. Run 'tauso setup-attract' to download it."
        )

    if not os.path.exists(expression_dir):
        raise FileNotFoundError(f"CRITICAL: Expression directory not found at {expression_dir}")

    # --- 2. MAP MATRIX ID -> GENE NAME ---
    logger.info("Loading metadata from %s to map Matrix IDs to Gene Names...", attract_csv_path)
    meta_df = pd.read_csv(attract_csv_path)

    # Filter metadata to only include the matrices currently in your pwm_db
    valid_matrix_ids = set(pwm_db.keys())
    filtered_meta = meta_df[meta_df["Matrix_id"].isin(valid_matrix_ids)]

    if filtered_meta.empty:
        raise ValueError("CRITICAL: No matching matrices found in ATtRACT CSV for the provided pwm_db keys.")

    # Get the unique list of REAL Gene names
    target_gene_names = filtered_meta["Gene_name"].unique().tolist()
    logger.info("Mapped %d matrices to %d unique gene names.", len(valid_matrix_ids), len(target_gene_names))

    # --- 3. LOAD EXPRESSION DATA ---
    unique_cell_ids = list(set(cell_lines))
    logger.info("Loading expression for %d cell lines...", len(unique_cell_ids))

    # Use existing function to load raw data
    transcriptomes = load_cell_line_gene_expression(
        depmap_ids=unique_cell_ids,
        valid_genes=target_gene_names,
        expression_dir=expression_dir,
    )

    if not transcriptomes:
        raise ValueError(
            f"CRITICAL: No expression data loaded! Checked {len(unique_cell_ids)} cell lines "
            f"against {len(target_gene_names)} genes. Check your paths and gene name formats."
        )

    # --- 4. CREATE EXPRESSION MATRIX ---
    logger.info("Building lookup matrix...")
    expr_list = []

    for ach_id, df_expr in transcriptomes.items():
        if df_expr.empty:
            continue

        # Keep only relevant columns
        temp_df = df_expr[["Gene", "expression_TPM"]].copy()

        # "Break the name" Logic (Safety Check)
        # Handles "GENE (ID)" format if present
        if not temp_df.empty and str(temp_df["Gene"].iloc[0]).endswith(")"):
            temp_df["Gene"] = temp_df["Gene"].astype(str).str.split(" ").str[0]

        temp_df[CELL_LINE_DEPMAP] = ach_id
        expr_list.append(temp_df)

    if not expr_list:
        raise ValueError("CRITICAL: Expression data was loaded but became empty during processing.")

    # Concatenate and Pivot
    full_expr_df = pd.concat(expr_list)

    # Pivot: Index=CellLine, Cols=GeneName, Values=TPM
    expression_matrix = full_expr_df.pivot(index=CELL_LINE_DEPMAP, columns="Gene", values="expression_TPM")

    if expression_matrix.empty:
        raise ValueError("CRITICAL: Final expression matrix is empty!")

    expression_matrix = np.log1p(expression_matrix)

    logger.info("Success. Expression matrix ready: %s", expression_matrix.shape)
    return expression_matrix
