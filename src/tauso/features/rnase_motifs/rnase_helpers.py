import json
from functools import lru_cache
from pathlib import Path

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


#################################################################################################
def score_window(subseq: str, weights_df) -> float:
    if not subseq:
        return 0.0

    score = 0.0
    for i, base in enumerate(subseq):
        if base not in weights_df.index:
            raise ValueError(f"Invalid base '{base}' at position {i}")
        col = f"Pos_{i + 1}"
        if col in weights_df.columns:
            score += weights_df.loc[base, col]
    return score / len(subseq)


def score_window_dict(subseq: str, weights: dict) -> float:
    if not subseq:
        return 0.0

    score = 0.0
    for i, base in enumerate(subseq):
        # Direct list access: O(1) and instant
        # weights['A'][0] corresponds to position 1
        try:
            score += weights[base][i]
        except (KeyError, IndexError):
            # Handle invalid base or sequence longer than weights
            pass

    return score / len(subseq)


#################################################################################################
def score_window_dinuc(subseq: str, weights_df) -> float:
    """
    Computes a dinucleotide-based score for a given subsequence using position-specific weights.
    Protects against empty subsequence and short lengths.
    """
    if not subseq or len(subseq) < 2:
        return 0.0

    score = 0.0
    for i in range(len(subseq) - 1):
        dimer = subseq[i : i + 2]
        col = f"Pos_{i + 1}"
        if dimer in weights_df.index and col in weights_df.columns:
            score += weights_df.loc[dimer, col]

    return score / len(subseq)


import numpy as np
from numba import njit

CHAR_TO_INT = {"A": 0, "C": 1, "G": 2, "T": 3}

# Cache to store converted matrices so we don't rebuild them on every single row
_WEIGHT_MATRIX_CACHE = {}


def get_numba_weights_matrix(weights_dict):
    """Converts the dictionary to a Numba-friendly NumPy matrix and caches it."""
    # Create a hashable key from the dict to use as a cache lookup
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
    """The ultra-fast math loop, stripped of strings and dicts."""
    L = len(seq_ints)
    if L < 2:
        return 0.0

    score = 0.0
    max_pos = weights_matrix.shape[1]

    for i in range(L - 1):
        # Math-based row lookup: (first_char * 4) + second_char
        dimer_idx = (seq_ints[i] * 4) + seq_ints[i + 1]

        if i < max_pos:
            score += weights_matrix[dimer_idx, i]

    return score / L


#######################################################################################################
def check_motif_presence(aso_sequence: str, motif: str) -> float:
    """
    Checks whether a specific motif is present in the ASO sequence.

    Parameters:
    -----------
    aso_sequence : str
        The DNA ASO sequence (A/C/G/T)

    motif : str
        A motif (sub-sequence) to search for (e.g., 'TCCC')

    Returns:
    --------
    float
        1.0 if the motif is found in the sequence, otherwise 0.0
        ---------
    good motif examples: 'TCCC', 'GGGA', 'CGCG', 'AGGA', 'TGCC' , 'CCCG'
    bad motif examples: 'TTTT', 'AAAA', 'CCCC', 'GGGG', 'TTAA', 'GCGC'
    """

    seq = aso_sequence.upper()
    motif = motif.upper()
    return 1.0 if motif in seq else 0.0


def compute_rnaseh1_score(aso_sequence: str, weights: dict, window_start: int = None) -> float:
    """
    Computes the RNase H1 cleavage score for a given ASO sequence using single-nucleotide weights.
    Refactored to use raw dictionary lookups for high performance.
    """

    # 1. Determine max_window from the weights dictionary
    # We assume all nucleotide lists in 'weights' are the same length.
    if not weights:
        return 0.0

    # Get the length of the weight vector from the first key (e.g., 'A')
    max_window = len(next(iter(weights.values())))

    # 2. Preprocess Sequence
    seq = aso_sequence.upper().replace("U", "T")
    L = len(seq)

    if L == 0:
        return 0.0

    # 3. Logic Branching

    # Case A: Sequence is shorter than the weight window -> Score full sequence
    if L < max_window:
        return score_window_dict(seq, weights)

    # Case B: Explicit window_start provided
    if window_start is not None:
        if window_start < 0 or window_start >= L:
            return 0.0
        # Slice from start. Python handles slicing past end automatically.
        subseq = seq[window_start : window_start + max_window]
        return score_window_dict(subseq, weights)

    # Case C: Default behavior (Center-aligned)
    if L % 2 == 0:
        # Even length: One central window
        start = (L - max_window) // 2
        subseq = seq[start : start + max_window]
        return score_window_dict(subseq, weights)
    else:
        # Odd length: Average of two central windows
        offset1 = (L - max_window) // 2
        offset2 = offset1 + 1

        s1 = score_window_dict(seq[offset1 : offset1 + max_window], weights)
        s2 = score_window_dict(seq[offset2 : offset2 + max_window], weights)
        return (s1 + s2) / 2


def compute_rnaseh1_dinucleotide_score(aso_sequence: str, dinuc_weights: dict, window_start: int = None) -> float:
    """
    Computes the RNase H1 cleavage score for a given ASO sequence using dinucleotide weights.
    (Optimized Numba Array Version)
    """
    if not dinuc_weights or not aso_sequence:
        return 0.0

    # 1. Get cached Numba matrix and max_window (instant lookup)
    weights_matrix, max_window = get_numba_weights_matrix(dinuc_weights)

    # 2. Preprocess Sequence into Integers ONCE
    seq_str = aso_sequence.upper().replace("U", "T")
    L = len(seq_str)

    if L == 0:
        return 0.0

    # Convert string to a NumPy array of ints. Unknown characters default to 0 ('A') to prevent crashing.
    seq_ints = np.array([CHAR_TO_INT.get(c, 0) for c in seq_str], dtype=np.int8)

    # 3. Logic Branching

    # Case A: Sequence is shorter than the weight window -> Score full sequence
    if L < max_window:
        return score_window_dinuc_dict(seq_ints, weights_matrix)

    # Case B: Explicit window_start provided
    if window_start is not None:
        if window_start < 0 or window_start >= L:
            return 0.0
        # Slice the numpy array (much faster than string slicing)
        subseq_ints = seq_ints[window_start : window_start + max_window]
        return score_window_dinuc_dict(subseq_ints, weights_matrix)

    # Case C: Default behavior (Center-aligned)
    if L % 2 == 0:
        start = (L - max_window) // 2
        subseq_ints = seq_ints[start : start + max_window]
        return score_window_dinuc_dict(subseq_ints, weights_matrix)
    else:
        offset1 = (L - max_window) // 2
        offset2 = offset1 + 1

        s1 = score_window_dinuc_dict(seq_ints[offset1 : offset1 + max_window], weights_matrix)
        s2 = score_window_dinuc_dict(seq_ints[offset2 : offset2 + max_window], weights_matrix)
        return (s1 + s2) / 2


# ==================== Global wrapper for RNaseH features ====================

# Default best window_start maps (length -> position)
BEST_WIN_NT = {
    "R4a": {
        10: 0,
        11: 0,
        12: 0,
        13: 0,
        14: 0,
        15: 0,
        16: 0,
        17: 2,
        18: 2,
        19: 4,
        20: 3,
        21: 0,
        22: 0,
        25: 0,
    },
    "R4b": {
        10: 0,
        11: 0,
        12: 0,
        13: 0,
        14: 0,
        15: 0,
        16: 0,
        17: 0,
        18: 0,
        19: 1,
        20: 3,
        21: 0,
        22: 0,
        25: 0,
    },
    "R7": {
        10: 0,
        11: 0,
        12: 0,
        13: 0,
        14: 0,
        15: 0,
        16: 3,
        17: 2,
        18: 4,
        19: 1,
        20: 4,
        21: 0,
        22: 0,
        25: 0,
    },
}

BEST_WIN_DINUC = {
    "R4a_dinuc": {
        10: 0,
        11: 0,
        12: 0,
        13: 0,
        14: 0,
        15: 0,
        16: 3,
        17: 3,
        18: 2,
        19: 4,
        20: 6,
        21: 0,
        22: 0,
        25: 0,
    },
    "R4b_dinuc": {
        10: 0,
        11: 0,
        12: 0,
        13: 0,
        14: 0,
        15: 0,
        16: 0,
        17: 0,
        18: 3,
        19: 1,
        20: 3,
        21: 0,
        22: 0,
        25: 0,
    },
    "R7_dinuc": {
        10: 0,
        11: 0,
        12: 0,
        13: 0,
        14: 0,
        15: 0,
        16: 0,
        17: 3,
        18: 4,
        19: 6,
        20: 3,
        21: 0,
        22: 0,
        25: 0,
    },
}

BEST_WIN_KREL_NT = {
    "R4a_krel": {
        10: 0,
        11: 0,
        12: 0,
        13: 0,
        14: 0,
        15: 0,
        16: 0,
        17: 3,
        18: 2,
        19: 4,
        20: 3,
        21: 0,
        22: 0,
        25: 0,
    },
    "R4b_krel": {
        10: 0,
        11: 0,
        12: 0,
        13: 0,
        14: 0,
        15: 0,
        16: 0,
        17: 1,
        18: 3,
        19: 1,
        20: 3,
        21: 0,
        22: 0,
        25: 0,
    },
    "R7_krel": {
        10: 0,
        11: 0,
        12: 0,
        13: 0,
        14: 0,
        15: 0,
        16: 3,
        17: 2,
        18: 4,
        19: 6,
        20: 4,
        21: 0,
        22: 0,
        25: 0,
    },
}

BEST_WIN_KREL_DINUC = {
    "R4a_krel_dinuc": {
        10: 0,
        11: 0,
        12: 0,
        13: 0,
        14: 0,
        15: 0,
        16: 2,
        17: 1,
        18: 2,
        19: 4,
        20: 6,
        21: 0,
        22: 0,
        25: 0,
    },
    "R4b_krel_dinuc": {
        10: 0,
        11: 0,
        12: 0,
        13: 0,
        14: 0,
        15: 0,
        16: 0,
        17: 0,
        18: 3,
        19: 1,
        20: 3,
        21: 0,
        22: 0,
        25: 0,
    },
    "R7_krel_dinuc": {
        10: 0,
        11: 0,
        12: 0,
        13: 0,
        14: 0,
        15: 0,
        16: 0,
        17: 3,
        18: 4,
        19: 6,
        20: 3,
        21: 0,
        22: 0,
        25: 0,
    },
}


def scan_constrained_window(target_seq: str, weights: dict, gap_start: int, gap_end: int) -> float:
    """
    Scans the target sequence with a sliding window defined by the weights.
    Enforces 'Best Effort' overlap constraints regarding the DNA gap.

    Args:
        target_seq: The full RNA target sequence (reverse complement of ASO).
        weights: Dictionary of nucleotide weights.
        gap_start: Start index of the DNA gap in the target sequence.
        gap_end: End index of the DNA gap in the target sequence.
    """
    window_size = len(next(iter(weights.values())))
    gap_len = gap_end - gap_start
    L = len(target_seq)

    valid_scores = []

    for i in range(L - window_size + 1):
        win_start = i
        win_end = i + window_size

        # Constraint Check:
        # 1. Big Window (Window > Gap): Window must fully engulf the gap.
        # 2. Small Window (Window <= Gap): Window must be fully inside the gap.
        if window_size > gap_len:
            is_valid = (win_start <= gap_start) and (win_end >= gap_end)
        else:
            is_valid = (win_start >= gap_start) and (win_end <= gap_end)

        if is_valid:
            score = compute_rnaseh1_score(target_seq, weights, window_start=i)
            valid_scores.append(score)

    return max(valid_scores) if valid_scores else 0.0
