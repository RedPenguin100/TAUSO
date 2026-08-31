"""Direct access to ViennaRNA's MFE dynamic-programming matrices.

`RNA.fold` answers "what is the MFE of this string". This module answers "what is the
MFE of positions a..b of this string", which ViennaRNA exposes no entry point for, by
folding the string once and reading every answer back out of the shared matrices.

That works because ``c[i][j]`` -- the best energy of the subsequence i..j given that i
pairs j -- is a function of that subsequence alone, so a single fold already contains
the interior of every window inside it. Only the exterior loop is context-dependent: a
stem sitting at a window edge has no dangling neighbour to stack on. `_exterior_loop_mfe`
recomputes just that part, per window, against the shared `c`.

Two things here sit lower than the rest of tauso. The `c` matrix is read through its raw
pointer, because copying it element-by-element out of the SWIG wrapper costs more than
the fold itself; and the exterior-loop recursion is a numba kernel, because it runs once
per window. `_check_matrix_layout` runs at import and fails loudly if ViennaRNA's memory
layout ever stops matching what is assumed here.
"""

import ctypes
from functools import lru_cache

import numpy as np
import ViennaRNA as RNA
from numba import njit

# ViennaRNA's "no valid structure" sentinel and its minimum hairpin loop size.
INF = 10000000
TURN = 3

# ViennaRNA hands its DP matrices to Python as `struct var_array { size_t length;
# T *data; unsigned int type; }`, so the payload pointer sits one size_t in.
_VAR_ARRAY_DATA_OFFSET = ctypes.sizeof(ctypes.c_size_t)

_BASE_CODE = np.zeros(128, dtype=np.int8)
for _base, _code in (("A", 1), ("C", 2), ("G", 3), ("U", 4)):
    _BASE_CODE[ord(_base)] = _code

# ViennaRNA's BP_pair: pair type of (5' base, 3' base), 0 where no pair is allowed.
_PAIR_TYPE = np.zeros((5, 5), dtype=np.int8)
for _five, _three, _pair_type in ((1, 4, 5), (2, 3, 1), (3, 2, 2), (3, 4, 3), (4, 1, 6), (4, 3, 4)):
    _PAIR_TYPE[_five, _three] = _pair_type


def _exterior_stem_energies():
    """E_ExtLoop for every (pair type, 5' neighbour, 3' neighbour); -1 means no neighbour."""
    params = RNA.param()
    table = np.full((8, 6, 6), INF, dtype=np.int32)
    for pair_type in range(1, 8):
        for five in range(-1, 5):
            for three in range(-1, 5):
                table[pair_type, five + 1, three + 1] = RNA.E_ExtLoop(pair_type, five, three, params)
    return table


_EXT_STEM = _exterior_stem_energies()


def _matrix_data(fold_compound, length):
    """Zero-copy view of ViennaRNA's `c` matrix, the base for the dynamic programming, read through its payload pointer.

    Indexing follows ViennaRNA's column-wise triangular layout,
    ``c[j * (j - 1) // 2 + i]`` for 1-based i < j.
    """
    address = ctypes.c_void_p.from_address(int(fold_compound.matrices.c.this) + _VAR_ARRAY_DATA_OFFSET).value
    return np.ctypeslib.as_array(
        ctypes.cast(address, ctypes.POINTER(ctypes.c_int32)),
        shape=(length * (length + 1) // 2 + 1,),
    )


# Small on purpose: each entry holds a full set of ViennaRNA DP matrices, not a float.
@lru_cache(maxsize=8)
def _folded(sequence, max_bp_span):
    """Fill ViennaRNA's MFE matrices for `sequence` and expose its `c` matrix.

    Returns (fold_compound, c, encoding). The fold_compound is carried along to keep the
    memory `c` views alive; `encoding` is 1-based to match `c`'s indexing.
    """
    model = RNA.md()
    model.max_bp_span = max_bp_span
    fold_compound = RNA.fold_compound(sequence, model)
    fold_compound.mfe()

    length = len(sequence)
    encoding = np.zeros(length + 1, dtype=np.int8)
    encoding[1:] = _BASE_CODE[np.frombuffer(sequence.encode(), dtype=np.uint8)]
    return fold_compound, _matrix_data(fold_compound, length), encoding


# TODO: add this function to ViennaRNA
# `_PAIR_TYPE` is read as a global but `_EXT_STEM` is passed in, on purpose. `cache=True`
# freezes a global's values into the on-disk cache, and the cache key is only this function's
# code, so a global that later changes is silently ignored. `_PAIR_TYPE` is a fixed lookup
# table; `_EXT_STEM` comes from `RNA.param()` and moves with the installed ViennaRNA.
@njit(cache=True)
def _exterior_loop_mfe(c, encoding, ext_stem, start, end):
    """MFE in dekacal of `encoding[start..end]` (1-based, inclusive), given a filled `c`."""
    width = end - start + 1
    best_upto = np.zeros(width + 1, dtype=np.int32)
    for offset in range(1, width + 1):
        j = start + offset - 1
        best = best_upto[offset - 1]
        column = (j * (j - 1)) // 2
        three_prime = -1 if j == end else encoding[j + 1]
        for i_offset in range(1, offset - TURN):
            i = start + i_offset - 1
            paired = c[column + i]
            if paired >= INF:
                continue
            pair_type = _PAIR_TYPE[encoding[i], encoding[j]]
            if pair_type == 0:
                continue
            five_prime = -1 if i == start else encoding[i - 1]
            energy = best_upto[i_offset - 1] + paired + ext_stem[pair_type, five_prime + 1, three_prime + 1]
            if energy < best:
                best = energy
        best_upto[offset] = best
    return best_upto[width]


def window_mfe(sequence: str, start: int, window_size: int, max_bp_span: int):
    """MFE of ``sequence[start:start + window_size]`` in kcal/mol.

    `sequence` is folded once and cached, so overlapping windows share the structure
    they have in common. Nothing inside a window can pair further than the window is
    wide, which is what caps the fold's base-pair span by default; callers reading
    several window sizes out of one fold pass the widest of them as `max_bp_span`, so
    they all land on the same cache entry. The float32 round matches how ViennaRNA
    converts its integer dekacal energies on the way out.
    """
    if window_size > max_bp_span:
        raise ValueError("[ViennaFold] window_mfe: window size is bigger than max_bp_span")
    _, c, encoding = _folded(sequence, max_bp_span)
    dekacal = _exterior_loop_mfe(c, encoding, _EXT_STEM, start + 1, start + window_size)
    return float(np.float32(dekacal / 100.0))


def _check_matrix_layout():
    """Fail at import if ViennaRNA's var_array layout is not the one `_folded` assumes."""
    sequence = "GGGGAAAACCCCAUAUAUGGCCUUAAGGCAUCGAUCG"
    fold_compound, c, _ = _folded(sequence, len(sequence))
    reference = fold_compound.matrices.c
    if any(int(c[i]) != reference[i] for i in range(0, len(reference), 7)):
        raise RuntimeError(
            "ViennaRNA's DP matrix memory layout differs from the one tauso reads. "
            "Check `struct var_array` in the installed ViennaRNA interface."
        )
    _folded.cache_clear()


_check_matrix_layout()
