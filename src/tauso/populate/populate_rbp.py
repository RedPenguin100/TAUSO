import logging
import math

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
def _log_unbound_numba_core(flat_seq, offsets, rows, weights, log_unbound):
    """Add each row's log-probability that this PWM leaves every site unoccupied.

    The sum over gapless placements of log(1 - o), where o = 1/(1 + 2^-s) is the placement's
    occupancy. Working in log space keeps the value stable when occupancies approach 1.

    flat_seq: every sequence end to end, as ints (A=0, C=1, G=2, U/T=3).
    offsets: where each sequence starts and ends in flat_seq.
    rows: the rows to score, those sharing one background.
    weights: (motif_len, 4) log2-odds of the PPM against that background.
    log_unbound: accumulated into, since an RBP's PWMs multiply.
    """
    motif_len = weights.shape[0]

    for row in rows:
        start = offsets[row]
        seq_len = offsets[row + 1] - start

        total = 0.0
        for i in range(seq_len - motif_len + 1):
            score = 0.0
            for pos in range(motif_len):
                score += weights[pos, flat_seq[start + i + pos]]
            total += math.log1p(-_occupancy_from_log2_odds(score))  # log(1 - o)

        log_unbound[row] += total


BASE_INDEX = {"A": 0, "C": 1, "G": 2, "U": 3, "T": 3}


def encode_sequences(sequences):
    """(flat, offsets): the sequences end to end as PWM column indices, and their boundaries.

    Encoding once is what makes the scan cheap: every PWM reads the same sequences, and
    converting a string costs more than scanning it.
    """
    flat, offsets = [], np.zeros(len(sequences) + 1, dtype=np.int64)
    for i, sequence in enumerate(sequences):
        text = "" if sequence is None or pd.isna(sequence) else str(sequence).upper()
        unknown = sorted(set(text) - BASE_INDEX.keys())
        if unknown:
            raise ValueError(f"Unknown base(s) {unknown} in sequence {text!r}; only A/C/G/U/T are allowed.")
        flat.extend(BASE_INDEX[base] for base in text)
        offsets[i + 1] = len(flat)
    return np.array(flat, dtype=np.int8), offsets


def log_odds_weights(pwm_matrix, background_probs):
    """A PWM's per-position log2-odds against the background, (motif_len, 4)."""
    pwm = pwm_matrix.values if hasattr(pwm_matrix, "values") else np.asarray(pwm_matrix)
    return np.log2((pwm.astype(np.float64) + 1e-9) / np.asarray(background_probs, dtype=np.float64))


def motif_log_unbound_numba(sequence, pwm_matrix, background_probs=None):
    """Σ log(1 - o) over a PWM's gapless placements (the log-probability it occupies no site).
    NaN/empty sequences contribute 0 (an unoccupied factor). See _log_unbound_numba_core."""
    if background_probs is None:
        background_probs = np.array([0.25, 0.25, 0.25, 0.25], dtype=np.float64)

    flat, offsets = encode_sequences([sequence])
    log_unbound = np.zeros(1)
    _log_unbound_numba_core(flat, offsets, np.arange(1), log_odds_weights(pwm_matrix, background_probs), log_unbound)
    return float(log_unbound[0])


def process_rbp(task, flat_seq, offsets, background_groups):
    """Worker function: processes ONE RBP for ALL sequences.

    The per-RBP score is the probability that the protein occupies at least one site in the
    window -- a noisy-OR, 1 - prod(1 - o), over the placements of all of its PWMs. Rows are
    grouped by background, so a PWM's weights are built once per background, not once per row.
    """
    matrices = task["matrices"]
    col_name = task["col_name"]

    log_unbound = np.zeros(offsets.shape[0] - 1, dtype=np.float64)
    for matrix in matrices:
        for background, rows in background_groups:
            _log_unbound_numba_core(flat_seq, offsets, rows, log_odds_weights(matrix, background), log_unbound)

    return col_name, 1.0 - np.exp(log_unbound)


def populate_rbp_affinity_features(df, rbp_map, pwm_db, gene_to_data, sequence_col="flank_sequence_50", n_jobs=32):
    """
    Calculates the raw Affinity for each RBP (regardless of Expression).
    """
    flank_param = sequence_col.split("_")[-1]
    df = df.loc[:, ~df.columns.duplicated()].copy()  # Use .copy() to avoid SettingWithCopy warnings later

    # --- 1. OPTIMIZATION: encode every sequence once, since all PWMs read the same ones ---
    flat_seq, offsets = encode_sequences(df[sequence_col].tolist())
    n_rows = len(df)

    # --- 1b. Backgrounds are a property of the gene, so count bases once per gene, and group
    # the rows that share one: a PWM's weights are then built per background, not per row.
    default_bg = np.array([0.25, 0.25, 0.25, 0.25], dtype=np.float32)
    genes = df[CANONICAL_GENE_NAME].to_numpy()
    bg_by_gene = {
        g: (get_background_probs(gene_to_data[g].full_mrna) if g in gene_to_data else default_bg)
        for g in pd.unique(genes)
    }
    background_probs_arr = np.array([bg_by_gene[g] for g in genes], dtype=np.float32)
    unique_bg, bg_index = np.unique(background_probs_arr, axis=0, return_inverse=True)
    background_groups = [(unique_bg[i], np.flatnonzero(bg_index == i)) for i in range(unique_bg.shape[0])]

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
    # A worker takes about a tenth of a second to start and is handed the PWM database, which
    # only pays off once the scan is long enough to hide it: under a few hundred rows one
    # process beats any number of them.
    workers = 1 if n_rows < 500 else n_jobs
    results = Parallel(n_jobs=workers)(
        delayed(process_rbp)(task, flat_seq, offsets, background_groups)
        for task in tqdm(target_tasks, desc="Computing RBPs", disable=workers == 1 and n_rows < 500)
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
