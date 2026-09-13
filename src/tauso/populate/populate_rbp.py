import logging
import math

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from numba import njit
from scipy.stats import entropy
from tqdm import tqdm

logger = logging.getLogger(__name__)

COMPLEXITY_ROW_CHUNK = 20000
"""Rows reduced at once when summing across RBPs."""

from tauso.algorithms.genomic_context_windows import flank_sequence_column
from tauso.data.consts import CANONICAL_GENE_NAME
from tauso.features.rbp.rbp_features import get_background_probs
from tauso.util import BASE_INDEX


@njit(fastmath=True)
def _occupancy_from_log2_odds(score):
    """Bound a per-window PWM log2-odds score to its [0,1] occupancy probability."""
    return 1.0 / (1.0 + 2.0 ** (-score))


@njit(fastmath=True)
def _log_unbound_numba_core(seq_indices, weights):
    """Log-probability that this PWM leaves every site unoccupied over the sequence.

    Returns the sum over gapless placements of log(1 - o), where o = 1/(1 + 2^-s) is the
    placement's occupancy. Working in log space keeps the value stable when occupancies
    approach 1.

    seq_indices: ints, one per nucleotide (A=0, C=1, G=2, U/T=3); -1 is not allowed.
    weights: (motif_len, 4) log2-odds against the background, columns A/C/G/U.
    """
    seq_len = len(seq_indices)
    motif_len = weights.shape[0]

    if seq_len < motif_len:
        return 0.0

    log_unbound = 0.0

    for i in range(seq_len - motif_len + 1):
        score = 0.0
        for pos in range(motif_len):
            score += weights[pos, seq_indices[i + pos]]
        log_unbound += math.log1p(-_occupancy_from_log2_odds(score))  # log(1 - o)

    return log_unbound


def _encode_sequence(sequence):
    """Map A/C/G/U/T to PWM columns; missing sequences have no sites."""
    if sequence is None or pd.isna(sequence):
        return np.empty(0, dtype=np.int8)
    seq_str = str(sequence)
    seq_indices = np.array([BASE_INDEX.get(base, -1) for base in seq_str.upper()], dtype=np.int8)
    if (seq_indices == -1).any():
        unknown = sorted(set(seq_str.upper()) - {"A", "C", "G", "U", "T"})
        raise ValueError(f"Unknown base(s) {unknown} in sequence {seq_str!r}; only A/C/G/U/T are allowed.")
    return seq_indices


def motif_log_unbound_numba(sequence, pwm_matrix, background_probs=None):
    """Σ log(1 - o) over a PWM's gapless placements (the log-probability it occupies no site).
    NaN/empty sequences contribute 0 (an unoccupied factor). See _log_unbound_numba_core."""
    seq_indices = _encode_sequence(sequence)
    if not len(seq_indices):
        return 0.0
    if background_probs is None:
        background_probs = np.full(4, 0.25)
    weights = np.log2((np.asarray(pwm_matrix, dtype=np.float64) + 1e-9) / background_probs)
    return _log_unbound_numba_core(seq_indices, weights)


@njit(fastmath=True)
def _log_unbound_batch(flat_seq, offsets, rows, weights, out):
    """Scan rows sharing a background without copying their sequences."""
    for row in rows:
        out[row] += _log_unbound_numba_core(flat_seq[offsets[row] : offsets[row + 1]], weights)


def process_rbp(task, flat_seq, offsets, background_groups):
    """Worker function: processes ONE RBP for ALL sequences.

    The per-RBP score is the probability that the protein occupies at least one site in the
    window -- a noisy-OR, 1 - prod(1 - o), over the placements of all of its PWMs.
    """
    matrices = task["matrices"]
    col_name = task["col_name"]
    n_rows = len(offsets) - 1

    log_unbound = np.zeros(n_rows, dtype=np.float64)
    for matrix in matrices:
        pwm = np.asarray(matrix, dtype=np.float64)
        for background, rows in background_groups:
            weights = np.log2((pwm + 1e-9) / background)
            _log_unbound_batch(flat_seq, offsets, rows, weights, log_unbound)

    return col_name, 1.0 - np.exp(log_unbound)


def populate_rbp_affinity_features(df, rbp_map, pwm_db, gene_to_data, flank_size, n_jobs=32):
    """One affinity column per RBP, as a frame indexed like `df`.

    Reads the two columns it needs and returns only what it computed. The step runs 27th of
    28, when `df` is at its widest, and copying it costs more than the scan does.
    """
    sequences = df[flank_sequence_column(flank_size)].fillna("").astype(str)
    genes = df[CANONICAL_GENE_NAME]
    n_rows = len(sequences)

    # Calculate backgrounds once per gene, then group rows sharing motif weights.
    default_bg = np.array([0.25, 0.25, 0.25, 0.25], dtype=np.float32)
    backgrounds = {
        g: get_background_probs(gene_to_data[g].full_mrna) if g in gene_to_data else default_bg for g in genes.unique()
    }
    background_probs_arr = np.array(genes.map(backgrounds).tolist(), dtype=np.float32).reshape(-1, 4)
    unique, inverse = np.unique(background_probs_arr, axis=0, return_inverse=True)
    background_groups = [(bg.astype(np.float64), np.flatnonzero(inverse == i)) for i, bg in enumerate(unique)]

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
                "col_name": f"rbp_{rbp.lower()}_aff_{flank_size}",
            }
        )

    if not target_tasks:
        logger.warning("No valid RBP tasks found in PWM DB.")
        return pd.DataFrame(index=sequences.index)

    logger.info("Calculating affinity features for %d RBPs on %d rows...", len(target_tasks), n_rows)

    # Encode once for all motifs; offsets delimit each row in the shared array.
    encoded = [_encode_sequence(sequence) for sequence in sequences]
    offsets = np.concatenate(([0], np.cumsum([len(sequence) for sequence in encoded], dtype=np.int64)))
    flat_seq = np.concatenate(encoded) if encoded else np.empty(0, dtype=np.int8)

    # --- 3. EXECUTION: Parallelize over RBPs, not Rows ---
    # Small batches cost less to scan than to dispatch to workers.
    results = Parallel(n_jobs=1 if n_rows < 500 else n_jobs)(
        delayed(process_rbp)(task, flat_seq, offsets, background_groups)
        for task in tqdm(target_tasks, desc="Computing RBPs")
    )

    logger.info("Done. Added %d affinity features.", len(results))
    return pd.DataFrame(dict(results), index=sequences.index)


def _total_and_diversity(df, feature_cols):
    """The per-row sum across `feature_cols`, and the entropy of the row once normalised.

    Both reduce along a row, so rows are taken a chunk at a time: the matrix and its
    normalised copy never exist for the whole frame at once.
    """
    total_scores = np.empty(len(df))
    diversity = np.empty(len(df))
    positions = [df.columns.get_loc(column) for column in feature_cols]

    for start in range(0, len(df), COMPLEXITY_ROW_CHUNK):
        stop = start + COMPLEXITY_ROW_CHUNK
        # Filling NaNs with 0 is crucial if some lookups failed.
        matrix = df.iloc[start:stop, positions].to_numpy(dtype=np.float64, na_value=0.0)

        chunk_total = matrix.sum(axis=1)
        total_scores[start:stop] = chunk_total

        # Normalize rows to sum to 1 to treat as probabilities.
        # Avoid division by zero for rows with 0 interaction.
        with np.errstate(divide="ignore", invalid="ignore"):
            matrix /= chunk_total[:, None]
            # Replace NaNs (from 0/0) with 0
            np.nan_to_num(matrix, copy=False)

        # Entropy, base e by default.
        diversity[start:stop] = entropy(matrix, axis=1)

    return total_scores, diversity


def populate_complexity_features(df, feature_cols, suffix, type="generic"):
    """
    Calculates role-independent global features:
    1. Total Interaction Load (Sum)
    2. Global Diversity (Entropy)
    """
    total_col = f"rbp_interaction_total_{suffix}_{type}"
    div_col = f"rbp_diversity_global_{suffix}_{type}"

    df[total_col], df[div_col] = _total_and_diversity(df, feature_cols)

    return df, [total_col, div_col]
