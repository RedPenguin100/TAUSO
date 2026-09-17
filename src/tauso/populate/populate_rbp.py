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
from tauso.data.consts import STRUCTURE_SENSE_START
from tauso.util import BASE_INDEX

BACKGROUND = 0.25
"""The nucleotide background the PWM log-odds are taken against. Uniform, so a motif scores
the same wherever it sits; the composition of the transcript around it plays no part."""


@njit(fastmath=True)
def _occupancy_from_log2_odds(score):
    """Bound a per-window PWM log2-odds score to its [0,1] occupancy probability."""
    return 1.0 / (1.0 + 2.0 ** (-score))


@njit(fastmath=True, nogil=True)
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


def _base_columns():
    """A byte's PWM column, -1 for anything that is not a base. Either case of a base reads the same."""
    table = np.full(256, -1, dtype=np.int8)
    for base, column in BASE_INDEX.items():
        table[ord(base)] = column
        table[ord(base.lower())] = column
    return table


_BASE_COLUMN = _base_columns()


def encode_sequences(sequences):
    """Glue the sequences into one array of PWM column indices.

    Returns (flat_seq, offsets): row i is flat_seq[offsets[i]:offsets[i + 1]].
    """
    sequences = list(sequences)
    offsets = np.cumsum([0, *map(len, sequences)], dtype=np.int64)
    letters = np.frombuffer("".join(sequences).encode("ascii"), dtype=np.uint8)
    flat_seq = _BASE_COLUMN[letters]
    if (flat_seq == -1).any():
        unknown = sorted({chr(b) for b in letters[flat_seq == -1]})
        raise ValueError(f"Unknown base(s) {unknown}; only A/C/G/U/T are allowed.")
    return flat_seq, offsets


@njit(fastmath=True, nogil=True)
def _log_unbound_batch(flat_seq, offsets, weights, out):
    """Scan every row without copying its sequence."""
    for row in range(len(out)):
        out[row] += _log_unbound_numba_core(flat_seq[offsets[row] : offsets[row + 1]], weights)


def process_rbp(task, flat_seq, offsets):
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
        weights = np.log2((pwm + 1e-9) / BACKGROUND)
        _log_unbound_batch(flat_seq, offsets, weights, log_unbound)

    return col_name, 1.0 - np.exp(log_unbound)


def populate_rbp_affinity_features(df, rbp_map, pwm_db, flank_size, n_jobs=32):
    """One affinity column per RBP, as a frame indexed like `df`.

    Reads the one column it needs and returns only what it computed. The step runs 27th of
    28, when `df` is at its widest, and copying it costs more than the scan does.
    """
    placed = df[STRUCTURE_SENSE_START] != -1
    sequences = df.loc[placed, flank_sequence_column(flank_size)]
    n_rows = len(sequences)

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
        return pd.DataFrame(index=df.index)

    logger.info("Calculating affinity features for %d RBPs on %d rows...", len(target_tasks), n_rows)

    # Encode once for all motifs; offsets delimit each row in the shared array.
    flat_seq, offsets = encode_sequences(sequences)

    # --- 3. EXECUTION: Parallelize over RBPs, not Rows ---
    # Threads, not processes: the kernel runs without the GIL, so they scale the same,
    # read flat_seq in place, and cost no interpreter of their own.
    # Small batches cost less to scan than to dispatch to workers.
    results = Parallel(n_jobs=1 if n_rows < 500 else n_jobs, prefer="threads")(
        delayed(process_rbp)(task, flat_seq, offsets) for task in tqdm(target_tasks, desc="Computing RBPs")
    )

    logger.info("Done. Added %d affinity features.", len(results))
    return pd.DataFrame(dict(results), index=sequences.index).reindex(df.index)


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
        matrix = df.iloc[start:stop, positions].to_numpy(dtype=np.float64)

        chunk_total = matrix.sum(axis=1)
        total_scores[start:stop] = chunk_total

        # Normalize rows to sum to 1 to treat as probabilities. A row with no window is NaN
        # in every column and stays NaN here.
        with np.errstate(divide="ignore", invalid="ignore"):
            matrix /= chunk_total[:, None]
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
