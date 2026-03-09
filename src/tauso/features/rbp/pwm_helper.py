import os

import numpy as np
import pandas as pd

from ...data.consts import CELL_LINE_DEPMAP
from ...data.data import get_data_dir
from ..codon_usage.find_cai_reference import load_cell_line_gene_expression


def calculate_total_affinity(sequence, pwm_matrix, background_probs=None, debug=False):
    # FORCE ARRAY
    if hasattr(pwm_matrix, 'values'):
        pwm_matrix = pwm_matrix.values
    else:
        pwm_matrix = np.array(pwm_matrix)

    if debug:
        print(f"      [AFFINITY-INTERNAL] Matrix Shape: {pwm_matrix.shape}")
        print(f"      [AFFINITY-INTERNAL] Seq Length: {len(str(sequence))}")
        print(f"      [AFFINITY-INTERNAL] Sequence Sample: {str(sequence)[:10]}...")

    if pd.isna(sequence) or len(sequence) == 0:
        if debug: print("      [AFFINITY-INTERNAL] ❌ Empty/NaN Sequence")
        return 0.0

    if background_probs is None:
        background_probs = np.array([0.25, 0.25, 0.25, 0.25])

    # MATH
    epsilon = 1e-9
    weights = np.log2((pwm_matrix + epsilon) / background_probs)

    base_map = {'A': 0, 'C': 1, 'G': 2, 'U': 3, 'T': 3}
    seq_indices = [base_map.get(base, -1) for base in str(sequence).upper()]

    seq_len = len(seq_indices)
    motif_len = len(pwm_matrix)

    if seq_len < motif_len:
        if debug: print(f"      [AFFINITY-INTERNAL] ❌ Seq too short ({seq_len} < {motif_len})")
        return 0.0

    total_score = 0.0

    # SCAN
    for i in range(seq_len - motif_len + 1):
        window = seq_indices[i: i + motif_len]
        if -1 in window: continue

        score = 0.0
        for pos, base_idx in enumerate(window):
            score += weights[pos, base_idx]

        if score > 0:
            total_score += score

    if debug:
        print(f"      [AFFINITY-INTERNAL] ✅ Final Score: {total_score}")

    return total_score


def build_rbp_expression_matrix(
        df: pd.DataFrame,
        pwm_db: dict,
        cell_line_col: str = CELL_LINE_DEPMAP,
        data_dir: str = None
) -> pd.DataFrame:
    """
    Builds a matrix of RNA expression (TPM) for the RBPs present in the pwm_db
    across the cell lines present in the provided dataframe.

    Args:
        df: DataFrame containing a column specifying cell lines (e.g., 'Cell_Line_Depmap').
        pwm_db: Dictionary of PWMs where keys are ATtRACT Matrix IDs.
        cell_line_col: Name of the column in df containing cell line IDs.
        data_dir: Path to data directory. If None, uses tauso.data.data.get_data_dir().

    Returns:
        pd.DataFrame: Pivot table (Index=CellLine, Columns=Gene, Values=TPM).
    """
    # --- 1. SETUP & FAIL LOUDLY ---
    if data_dir is None:
        data_dir = get_data_dir()

    attract_csv_path = os.path.join(data_dir, "RBS_motifs_Homo_sapiens.csv")
    expression_dir = os.path.join(data_dir, "processed_expression")

    if not os.path.exists(attract_csv_path):
        raise FileNotFoundError(f"CRITICAL: ATtRACT metadata not found at {attract_csv_path}")

    if not os.path.exists(expression_dir):
        raise FileNotFoundError(f"CRITICAL: Expression directory not found at {expression_dir}")

    # --- 2. MAP MATRIX ID -> GENE NAME ---
    print(f"Loading metadata from {attract_csv_path} to map Matrix IDs to Gene Names...")
    meta_df = pd.read_csv(attract_csv_path)

    # Filter metadata to only include the matrices currently in your pwm_db
    valid_matrix_ids = set(pwm_db.keys())
    filtered_meta = meta_df[meta_df['Matrix_id'].isin(valid_matrix_ids)]

    if filtered_meta.empty:
        raise ValueError("CRITICAL: No matching matrices found in ATtRACT CSV for the provided pwm_db keys.")

    # Get the unique list of REAL Gene names
    target_gene_names = filtered_meta['Gene_name'].unique().tolist()
    print(f"Mapped {len(valid_matrix_ids)} matrices to {len(target_gene_names)} unique gene names.")

    # --- 3. LOAD EXPRESSION DATA ---
    unique_cell_ids = df[cell_line_col].unique().tolist()
    print(f"Loading expression for {len(unique_cell_ids)} cell lines...")

    # Use existing function to load raw data
    transcriptomes = load_cell_line_gene_expression(
        depmap_ids=unique_cell_ids,
        valid_genes=target_gene_names,
        expression_dir=expression_dir
    )

    if not transcriptomes:
        raise ValueError(
            f"CRITICAL: No expression data loaded! Checked {len(unique_cell_ids)} cell lines "
            f"against {len(target_gene_names)} genes. Check your paths and gene name formats."
        )

    # --- 4. CREATE EXPRESSION MATRIX ---
    print("Building lookup matrix...")
    expr_list = []

    for ach_id, df_expr in transcriptomes.items():
        if df_expr.empty:
            continue

        # Keep only relevant columns
        temp_df = df_expr[['Gene', 'expression_TPM']].copy()

        # "Break the name" Logic (Safety Check)
        # Handles "GENE (ID)" format if present
        if not temp_df.empty and str(temp_df['Gene'].iloc[0]).endswith(")"):
            temp_df['Gene'] = temp_df['Gene'].astype(str).str.split(" ").str[0]

        temp_df[cell_line_col] = ach_id
        expr_list.append(temp_df)

    if not expr_list:
        raise ValueError("CRITICAL: Expression data was loaded but became empty during processing.")

    # Concatenate and Pivot
    full_expr_df = pd.concat(expr_list)

    # Pivot: Index=CellLine, Cols=GeneName, Values=TPM
    expression_matrix = full_expr_df.pivot(
        index=cell_line_col,
        columns='Gene',
        values='expression_TPM'
    )

    if expression_matrix.empty:
        raise ValueError("CRITICAL: Final expression matrix is empty!")

    print(f"Success. Expression matrix ready: {expression_matrix.shape}")
    return expression_matrix
