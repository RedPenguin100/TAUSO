import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from numba import njit

_WEIGHTS_DIR = Path(__file__).resolve().parent / "weights"


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


def score_window_dict(subseq: str, weights: dict) -> float:
    """Mean per-position single-nucleotide weight for ``subseq``."""
    if not subseq:
        return 0.0

    score = 0.0
    for i, base in enumerate(subseq):
        try:
            score += weights[base][i]
        except (KeyError, IndexError):
            pass

    return score / len(subseq)


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
def score_window_dinuc_dict(seq_ints, weights_matrix):
    """Numba kernel: sum dimer weights along ``seq_ints``, normalized by length."""
    L = len(seq_ints)
    if L < 2:
        return 0.0

    score = 0.0
    max_pos = weights_matrix.shape[1]

    for i in range(L - 1):
        dimer_idx = (seq_ints[i] * 4) + seq_ints[i + 1]
        if i < max_pos:
            score += weights_matrix[dimer_idx, i]

    return score / L


def compute_rnaseh1_score(aso_sequence: str, weights: dict, window_start: int) -> float:
    """Score a single ``max_window``-wide slice of ``aso_sequence`` starting at ``window_start``.

    When the sequence is shorter than ``max_window`` the full sequence is scored
    (``window_start`` is ignored) — this short-seq fallback is preserved from the
    original implementation that produced the regression baseline.
    """
    if not weights or not aso_sequence:
        return 0.0

    max_window = len(next(iter(weights.values())))
    seq = aso_sequence.upper().replace("U", "T")
    L = len(seq)

    if L < max_window:
        return score_window_dict(seq, weights)

    if window_start < 0 or window_start >= L:
        return 0.0

    return score_window_dict(seq[window_start : window_start + max_window], weights)


def compute_rnaseh1_dinucleotide_score(aso_sequence: str, dinuc_weights: dict, window_start: int) -> float:
    """Dinucleotide counterpart of :func:`compute_rnaseh1_score`; uses the Numba kernel.

    Same short-sequence fallback as the single-nt scorer.
    """
    if not dinuc_weights or not aso_sequence:
        return 0.0

    weights_matrix, max_window = get_numba_weights_matrix(dinuc_weights)
    seq_str = aso_sequence.upper().replace("U", "T")
    L = len(seq_str)
    seq_ints = np.array([CHAR_TO_INT.get(c, 0) for c in seq_str], dtype=np.int8)

    if L < max_window:
        return score_window_dinuc_dict(seq_ints, weights_matrix)

    if window_start < 0 or window_start >= L:
        return 0.0

    return score_window_dinuc_dict(seq_ints[window_start : window_start + max_window], weights_matrix)


def scan_constrained_window(target_seq: str, weights: dict, gap_start: int, gap_end: int) -> float:
    """Max single-nt score over windows that satisfy a gap-overlap constraint.

    If the window is wider than the DNA gap, it must fully engulf the gap;
    otherwise the window must sit fully inside the gap.
    """
    window_size = len(next(iter(weights.values())))
    gap_len = gap_end - gap_start
    L = len(target_seq)

    valid_scores = []
    for i in range(L - window_size + 1):
        win_start = i
        win_end = i + window_size

        if window_size > gap_len:
            is_valid = (win_start <= gap_start) and (win_end >= gap_end)
        else:
            is_valid = (win_start >= gap_start) and (win_end <= gap_end)

        if is_valid:
            valid_scores.append(compute_rnaseh1_score(target_seq, weights, window_start=win_start))

    return max(valid_scores) if valid_scores else 0.0
