import logging

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from numba import njit
from scipy.stats import entropy
from tqdm import tqdm

logger = logging.getLogger(__name__)

from tauso.data.consts import CANONICAL_GENE_NAME
from tauso.features.rbp.rbp_features import get_background_probs


@njit(fastmath=True)
def _occupancy_from_log2_odds(score):
    """Bound a per-window PWM log2-odds score to its [0,1] occupancy probability."""
    return 1.0 / (1.0 + 2.0 ** (-score))


@njit(fastmath=True)
def _calculate_affinity_numba_core(seq_indices, pwm_matrix, background_probs):
    """Return the total motif-occupancy affinity of one PWM over one sequence.

    seq_indices: ints, one per nucleotide (A=0, C=1, G=2, U/T=3); -1 is not allowed.
    pwm_matrix: (motif_len, 4) PPM, columns P(A),P(C),P(G),P(U), each row summing to ~1.
    background_probs: null base frequencies [pA, pC, pG, pU].
    """
    seq_len = len(seq_indices)
    motif_len = pwm_matrix.shape[0]

    if seq_len < motif_len:
        return 0.0

    weights = np.log2((pwm_matrix + 1e-9) / background_probs)

    total_score = 0.0

    for i in range(seq_len - motif_len + 1):
        score = 0.0
        for pos in range(motif_len):
            score += weights[pos, seq_indices[i + pos]]
        total_score += _occupancy_from_log2_odds(score)

    return total_score


def calculate_total_affinity_numba(sequence, pwm_matrix, background_probs=None, debug=False):
    # 1. Type guard the matrix
    if hasattr(pwm_matrix, "values"):
        pwm_matrix = pwm_matrix.values
    else:
        pwm_matrix = np.asarray(pwm_matrix)

    # 2. Handle empty/NaN sequences
    if pd.isna(sequence) or sequence is None:
        return 0.0

    seq_str = str(sequence)
    if len(seq_str) == 0:
        return 0.0

    # 3. Handle background probs
    if background_probs is None:
        background_probs = np.array([0.25, 0.25, 0.25, 0.25], dtype=np.float64)
    else:
        background_probs = np.asarray(background_probs, dtype=np.float64)

    # 4. Map string to integers (Pure Python is fine here before handing off to Numba)
    base_map = {"A": 0, "C": 1, "G": 2, "U": 3, "T": 3}
    seq_indices = np.array([base_map.get(base, -1) for base in seq_str.upper()], dtype=np.int8)

    # Unknown bases are not allowed: the scorer needs the full sequence, so fail loudly.
    if (seq_indices == -1).any():
        unknown = sorted(set(seq_str.upper()) - {"A", "C", "G", "U", "T"})
        raise ValueError(f"Unknown base(s) {unknown} in sequence {seq_str!r}; only A/C/G/U/T are allowed.")

    # 5. Execute the compiled loop
    # We enforce float64 to ensure math precision matches your original numpy logic perfectly
    return _calculate_affinity_numba_core(seq_indices, pwm_matrix.astype(np.float64), background_probs)


def process_rbp(task, sequences, background_probs_arr):
    """
    Worker function: Processes ONE RBP for ALL sequences.
    Sums per-PWM affinities across every PWM that ATtRACT lists for this RBP.
    """
    matrices = task["matrices"]
    col_name = task["col_name"]
    n_rows = len(sequences)

    scores = np.zeros(n_rows, dtype=np.float64)
    for matrix in matrices:
        for i in range(n_rows):
            scores[i] += calculate_total_affinity_numba(sequences[i], matrix, background_probs_arr[i])

    return col_name, scores


def populate_rbp_affinity_features(df, rbp_map, pwm_db, gene_to_data, sequence_col="flank_sequence_50", n_jobs=32):
    """
    Calculates the raw Affinity for each RBP (regardless of Expression).
    """
    flank_param = sequence_col.split("_")[-1]
    df = df.loc[:, ~df.columns.duplicated()].copy()  # Use .copy() to avoid SettingWithCopy warnings later

    # Keep sequences as a list for now, as strings don't map well in NumPy
    sequences = df[sequence_col].fillna("").astype(str).tolist()
    n_rows = len(sequences)

    # --- 1. OPTIMIZATION: Create a 2D NumPy array for Background Probs ---
    # This allows Joblib to share memory across workers instead of pickling 180k tiny lists
    default_bg = np.array([0.25, 0.25, 0.25, 0.25], dtype=np.float32)
    bg_list = [
        get_background_probs(gene_to_data[g].full_mrna) if g in gene_to_data else default_bg
        for g in df[CANONICAL_GENE_NAME]
    ]
    background_probs_arr = np.array(bg_list, dtype=np.float32)

    # --- 2. FILTER & PREPARE RBP METADATA ---
    target_tasks = []
    for rbp, matrix_ids in rbp_map.items():
        valid_mids = [m for m in matrix_ids if m in pwm_db]
        if not valid_mids:
            continue

        target_tasks.append(
            {
                "name": rbp,
                "matrices": [pwm_db[m] for m in valid_mids],
                "col_name": f"rbp_{rbp.lower()}_aff_{flank_param}",
            }
        )

    if not target_tasks:
        logger.warning("No valid RBP tasks found in PWM DB.")
        return df, []

    logger.info("Calculating affinity features for %d RBPs on %d rows...", len(target_tasks), n_rows)

    # --- 3. EXECUTION: Parallelize over RBPs, not Rows ---
    # joblib uses 'loky' backend by default, which excels at memory mapping large arrays (background_probs_arr)
    results = Parallel(n_jobs=n_jobs)(
        delayed(process_rbp)(task, sequences, background_probs_arr)
        for task in tqdm(target_tasks, desc="Computing RBPs")
    )

    # --- 4. AGGREGATION & ASSIGNMENT ---
    logger.debug("Assigning columns...")

    # Bundle all new columns into a dictionary, then concat once
    new_cols_dict = {col_name: scores for col_name, scores in results}
    df = pd.concat([df, pd.DataFrame(new_cols_dict, index=df.index)], axis=1)
    new_col_names = list(new_cols_dict.keys())

    logger.info("Done. Added %d affinity features.", len(new_col_names))
    return df, new_col_names


def populate_complexity_features(df, feature_cols, suffix, type="generic"):
    """
    Calculates role-independent global features:
    1. Total Interaction Load (Sum)
    2. Global Diversity (Entropy)
    """
    # Create a matrix of the relevant columns (add epsilon to avoid div/0)
    # Filling NaNs with 0 is crucial if some lookups failed
    matrix = df[feature_cols].fillna(0.0).values

    # 1. Total Interaction (Sum of all RBPs)
    total_col = f"rbp_interaction_total_{suffix}_{type}"
    total_scores = matrix.sum(axis=1)
    df[total_col] = total_scores

    # 2. Global Diversity (Shannon Entropy)
    # Normalize rows to sum to 1 to treat as probabilities
    # Avoid division by zero for rows with 0 interaction
    with np.errstate(divide="ignore", invalid="ignore"):
        probs = matrix / total_scores[:, None]
        # Replace NaNs (from 0/0) with 0
        probs = np.nan_to_num(probs)

    div_col = f"rbp_diversity_global_{suffix}_{type}"
    # Calculate entropy (base e by default)
    df[div_col] = entropy(probs, axis=1)

    return df, [total_col, div_col]
