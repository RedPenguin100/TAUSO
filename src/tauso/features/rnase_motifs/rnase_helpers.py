import json
import logging
import math
from functools import lru_cache
from pathlib import Path

import numpy as np
from numba import njit

logger = logging.getLogger(__name__)

_WEIGHTS_DIR = Path(__file__).resolve().parent / "weights"


def _gap_has_unknown_bases(seq: str, gap_start: int, gap_end: int) -> bool:
    gap = seq[gap_start:gap_end].upper().replace("U", "T")
    return any(b not in "ACGT" for b in gap)


@lru_cache(maxsize=None)
def rnaseh1_dict(label: str):
    """Position-specific RNase H1 weights for the named experiment.

    Bundled JSON under ``weights/`` — one file per (experiment, encoding):
    ``R4a.json``, ``R4a_dinuc.json``, ``R7_krel_dinuc.json``, etc. Source:
    Kiełpiński et al. 2017 (NAR; PMC5728404), experiments R4a / R4b / R7.
    """
    path = _WEIGHTS_DIR / f"{label}.json"
    if not path.exists():
        raise ValueError(f"Unknown RNase H1 weights label: {label!r}")
    return json.loads(path.read_text())


CHAR_TO_INT = {"A": 0, "C": 1, "G": 2, "T": 3}

_WEIGHT_MATRIX_CACHE: dict = {}


def get_numba_weights_matrix(weights_dict):
    """Convert the dinuc weights dict to a (16, max_w) NumPy matrix; cached by content."""
    cache_key = tuple(sorted((k, tuple(v)) for k, v in weights_dict.items()))

    if cache_key not in _WEIGHT_MATRIX_CACHE:
        max_w = len(next(iter(weights_dict.values()))) if weights_dict else 0
        matrix = np.zeros((16, max_w), dtype=np.float64)

        for dimer, vals in weights_dict.items():
            if len(dimer) == 2 and dimer[0] in CHAR_TO_INT and dimer[1] in CHAR_TO_INT:
                row_idx = (CHAR_TO_INT[dimer[0]] * 4) + CHAR_TO_INT[dimer[1]]
                matrix[row_idx, : len(vals)] = vals

        _WEIGHT_MATRIX_CACHE[cache_key] = (matrix, max_w)

    return _WEIGHT_MATRIX_CACHE[cache_key]


@njit(fastmath=True)
def _score_dinuc_window_masked(seq_ints, weights_matrix, window_start, gap_start, gap_end):
    """Numba kernel: dinuc score over one window, with flank-bridging dimers zeroed."""
    window_size = weights_matrix.shape[1]
    L = len(seq_ints)
    score = 0.0
    for i in range(window_size - 1):
        pos1 = window_start + i
        pos2 = pos1 + 1
        if pos1 < gap_start or pos2 >= gap_end or pos2 >= L:
            continue
        dimer_idx = (seq_ints[pos1] * 4) + seq_ints[pos2]
        score += weights_matrix[dimer_idx, i]
    return score / window_size


def scan_constrained_window(target_seq: str, weights: dict, gap_start: int, gap_end: int) -> float:
    """Max single-nt PSSM score over windows overlapping the DNA gap.

    Positions outside ``[gap_start, gap_end)`` contribute 0 (sugar-modified flanks
    are not cleaved). Returns NaN with a warning if the gap has any non-ACGT base.
    """
    if _gap_has_unknown_bases(target_seq, gap_start, gap_end):
        logger.warning("Non-ACGT base in DNA gap [%d, %d) of %r; returning NaN", gap_start, gap_end, target_seq)
        return math.nan

    window_size = len(next(iter(weights.values())))
    gap_len = gap_end - gap_start
    L = len(target_seq)
    seq = target_seq.upper().replace("U", "T")

    valid_scores = []
    for i in range(L - window_size + 1):
        win_start = i
        win_end = i + window_size

        if window_size > gap_len:
            is_valid = (win_start <= gap_start) and (win_end >= gap_end)
        else:
            is_valid = (win_start >= gap_start) and (win_end <= gap_end)

        if not is_valid:
            continue

        score = 0.0
        for k in range(window_size):
            pos = win_start + k
            if pos < gap_start or pos >= gap_end:
                continue
            try:
                score += weights[seq[pos]][k]
            except (KeyError, IndexError):
                pass
        valid_scores.append(score / window_size)

    return max(valid_scores) if valid_scores else 0.0


def scan_constrained_window_dinuc(target_seq: str, dinuc_weights: dict, gap_start: int, gap_end: int) -> float:
    """Dinucleotide counterpart of :func:`scan_constrained_window`. Same overlap rule
    and flank-zeroing. Returns NaN with a warning if the gap has any non-ACGT base.
    """
    if _gap_has_unknown_bases(target_seq, gap_start, gap_end):
        logger.warning("Non-ACGT base in DNA gap [%d, %d) of %r; returning NaN", gap_start, gap_end, target_seq)
        return math.nan

    weights_matrix, window_size = get_numba_weights_matrix(dinuc_weights)
    gap_len = gap_end - gap_start
    L = len(target_seq)
    seq_str = target_seq.upper().replace("U", "T")
    seq_ints = np.array([CHAR_TO_INT.get(c, 0) for c in seq_str], dtype=np.int8)

    valid_scores = []
    for i in range(L - window_size + 1):
        win_start = i
        win_end = i + window_size

        if window_size > gap_len:
            is_valid = (win_start <= gap_start) and (win_end >= gap_end)
        else:
            is_valid = (win_start >= gap_start) and (win_end <= gap_end)

        if is_valid:
            valid_scores.append(_score_dinuc_window_masked(seq_ints, weights_matrix, win_start, gap_start, gap_end))

    return max(valid_scores) if valid_scores else 0.0
