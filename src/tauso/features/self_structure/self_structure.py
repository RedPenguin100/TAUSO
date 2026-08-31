"""What an ASO does to itself: the best hairpin it can fold into and the best homodimer it
can form with a second copy.

Both compete with binding the target, and both are properties of the oligo alone -- sequence
and chemistry, no target involved. Scored with the nearest-neighbour weights in `nn_tables`,
so 2'-MOE and cEt positions carry their own stacking rather than being read as deoxy.

The search is exhaustive over stems rather than a partition function: the reported value is
the single best structure, not an ensemble average.
"""

import numpy as np
from numba import njit

from ...util import DNA_BASES
from .nn_tables import (
    CET,
    CLASS_MAP,
    COMPLEMENT,
    DEOXY,
    HAIRPIN_LOOP,
    INIT_AT,
    INIT_GC,
    MAX_LOOP,
    MIN_LOOP,
    MOE,
    PARAMETERS,
)

FEATURE_NAMES = [
    "hp_dg",
    "hp_stem",
    "hp_modpairs",
    "dim_dg",
    "dim_helix",
    "dim_modpairs",
    "dim_minus_hp",
    "hp_dg_per_nt",
    "dim_dg_per_nt",
]

_A, _T = DNA_BASES.index("A"), DNA_BASES.index("T")


@njit(cache=True)
def _pairs(a, b):
    return COMPLEMENT[a] == b


@njit(cache=True)
def _complement(a):
    return COMPLEMENT[a]


@njit(cache=True)
def _initiation(a):
    """SantaLucia helix initiation: an A-T terminus frays more easily than G-C."""
    return INIT_AT if (a == _A or a == _T) else INIT_GC


@njit(cache=True)
def _stack(s1, s2, c1, c2, d1, d2, parameters, class_map):
    """dG of one stack; (s1, s2) is the 5'->3' dinucleotide, c and d the sugars of both strands."""
    cls = class_map[c1, c2, d1, d2, 0]
    if class_map[c1, c2, d1, d2, 1] == 1:
        return parameters[cls, _complement(s2), _complement(s1)]
    return parameters[cls, s1, s2]


@njit(cache=True)
def hairpin(seq, length, sugars, parameters, class_map):
    """Best stem-loop the oligo folds into. Returns (dG, stem length, modified positions paired)."""
    best, best_stem, best_modified = 0.0, 0, 0
    for i in range(length):
        for j in range(i + MIN_LOOP + 1, length):
            if not _pairs(seq[i], seq[j]):
                continue
            loop = j - i - 1
            if loop < MIN_LOOP or loop > MAX_LOOP:
                continue
            base = HAIRPIN_LOOP[loop]
            stack = 0.0
            modified = (1 if sugars[i] != DEOXY else 0) + (1 if sugars[j] != DEOXY else 0)
            for k in range(1, length):
                left, right = i - k, j + k
                if left < 0 or right >= length or not _pairs(seq[left], seq[right]):
                    break
                stack += _stack(
                    seq[left],
                    seq[left + 1],
                    sugars[left],
                    sugars[left + 1],
                    sugars[right],
                    sugars[right - 1],
                    parameters,
                    class_map,
                )
                modified += (1 if sugars[left] != DEOXY else 0) + (1 if sugars[right] != DEOXY else 0)
                energy = base + stack + _initiation(seq[left])
                if energy < best:
                    best, best_stem, best_modified = energy, k + 1, modified
    return best, best_stem, best_modified


@njit(cache=True)
def dimer(seq, length, sugars, parameters, class_map):
    """Best duplex with a second copy of itself, over every register.

    Returns (dG, helix length, modified positions paired).
    """
    best, best_helix, best_modified = 0.0, 0, 0
    for shift in range(-(length - 1), length):
        i = 0
        while i < length:
            j = length - 1 - i - shift
            if 0 <= j < length and _pairs(seq[i], seq[j]):
                start = i
                while i + 1 < length:
                    following = length - 2 - i - shift
                    if 0 <= following < length and _pairs(seq[i + 1], seq[following]):
                        i += 1
                    else:
                        break
                end = i
                if end > start:
                    energy = _initiation(seq[start]) + _initiation(seq[end])
                    modified = 0
                    for t in range(start, end):
                        partner = length - 1 - t - shift
                        energy += _stack(
                            seq[t],
                            seq[t + 1],
                            sugars[t],
                            sugars[t + 1],
                            sugars[partner],
                            sugars[partner - 1],
                            parameters,
                            class_map,
                        )
                    for t in range(start, end + 1):
                        partner = length - 1 - t - shift
                        modified += (1 if sugars[t] != DEOXY else 0) + (1 if sugars[partner] != DEOXY else 0)
                    if energy < best:
                        best, best_helix, best_modified = energy, end - start + 1, modified
            i += 1
    return best, best_helix, best_modified


@njit(cache=True)
def _score_all(seqs, lengths, sugars, parameters, class_map):
    out = np.zeros((seqs.shape[0], 6), dtype=np.float64)
    for r in range(seqs.shape[0]):
        length = lengths[r]
        out[r, 0], out[r, 1], out[r, 2] = hairpin(seqs[r], length, sugars[r], parameters, class_map)
        out[r, 3], out[r, 4], out[r, 5] = dimer(seqs[r], length, sugars[r], parameters, class_map)
    return out


def encode(sequences, chemical_patterns):
    """Sequences and sugars as integer arrays, sized to the longest oligo given.

    A `chemical_pattern` letter of M is 2'-MOE and C is cEt; anything else is deoxy. Oligos
    carrying a letter outside ACGT are left at length zero and come back unscored, since a
    base the weights do not cover cannot be stacked.
    """
    index = {base: i for i, base in enumerate(DNA_BASES)}
    count = len(sequences)
    width = max((len(str(s)) for s in sequences), default=0)
    seqs = np.zeros((count, width), dtype=np.int8)
    sugars = np.zeros((count, width), dtype=np.int8)
    lengths = np.zeros(count, dtype=np.int64)
    for r, (sequence, pattern) in enumerate(zip(sequences, chemical_patterns)):
        sequence = str(sequence).upper()
        if any(base not in index for base in sequence):
            continue
        lengths[r] = len(sequence)
        for k, base in enumerate(sequence):
            seqs[r, k] = index[base]
        pattern = str(pattern)
        for k in range(min(len(pattern), len(sequence))):
            sugars[r, k] = MOE if pattern[k] == "M" else (CET if pattern[k] == "C" else DEOXY)
    return seqs, lengths, sugars


def calculate_self_structure(sequences, chemical_patterns):
    """The nine self-structure features, one row per oligo, in `FEATURE_NAMES` order.

    Unscorable oligos come back NaN rather than zero, so they are not read as a free oligo.
    """
    seqs, lengths, sugars = encode(sequences, chemical_patterns)
    raw = _score_all(seqs, lengths, sugars, PARAMETERS, CLASS_MAP)
    scored = lengths > 0
    per_nt = np.where(scored, np.maximum(lengths, 1), np.nan)
    out = {
        "hp_dg": raw[:, 0],
        "hp_stem": raw[:, 1],
        "hp_modpairs": raw[:, 2],
        "dim_dg": raw[:, 3],
        "dim_helix": raw[:, 4],
        "dim_modpairs": raw[:, 5],
        "dim_minus_hp": raw[:, 3] - raw[:, 0],
        "hp_dg_per_nt": raw[:, 0] / per_nt,
        "dim_dg_per_nt": raw[:, 3] / per_nt,
    }
    for name in FEATURE_NAMES:
        out[name] = np.where(scored, out[name], np.nan)
    return out
