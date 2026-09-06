"""What an oligo does with itself: folds into a hairpin, or pairs with a second copy.

Three searches, all over the oligo alone -- no target, no locus, no fold of anything else.

    hairpin        the best stem-loop the oligo folds into
    dimer          the best Watson-Crick helix between the oligo and a copy of it, over
                   every register
    whole overlap  the energy of the entire overlapping region at each register, with
                   mismatched pairs paying their fitted cost instead of ending the helix

The first two ask "what is the best matched stretch"; the third asks "what does the whole
duplex cost", which is a different question and needs mismatch weights to be answerable at
all. A typical oligo's best register is only about 72% Watson-Crick, so the third search
reaches structures the first two cannot describe.

Mismatches are crossed only in the whole-overlap search. Allowing them in the hairpin and
dimer searches was measured to make those features worse: nearly every oligo then finds some
long partially-paired stretch, and the feature stops separating oligos that genuinely
self-pair from ones that do not.
"""

import numpy as np
from numba import njit

from ...util import DNA_BASES, normalize_dna
from .weights import (
    COMPLEMENT,
    DEOXY,
    DINUCLEOTIDE_ENERGY,
    MISMATCH_ENERGY,
    TANDEM_MISMATCH_ENERGY,
)
from .santalucia import (
    HAIRPIN_LOOP_PENALTY,
    INITIATION_AT,
    INITIATION_GC,
    MAX_LOOP,
    MIN_LOOP,
)

SUGARS = {"M": 1, "C": 2}
"""`chemical_pattern` letters. Deoxy is anything else, which is how the column writes it."""

_A, _T = DNA_BASES.index("A"), DNA_BASES.index("T")

HAIRPIN_NAMES = ["hp_dg", "hp_stem", "hp_modpairs", "hp_dg_per_nt"]
PERFECT_NAMES = [
    "homodimer_perfect_dg",
    "homodimer_perfect_helix",
    "homodimer_perfect_modpairs",
    "homodimer_perfect_dg_per_nt",
    "homodimer_perfect_minus_hp",
]
MISMATCHED_NAMES = [
    "homodimer_mismatched_dg",
    "homodimer_mismatched_len",
    "homodimer_mismatched_nmm",
    "homodimer_mismatched_dg_per_nt",
    "homodimer_mismatched_full_dg",
    "homodimer_mismatched_full_nmm",
]
FEATURE_NAMES = HAIRPIN_NAMES + PERFECT_NAMES + MISMATCHED_NAMES


@njit(cache=True)
def _is_paired(first, second):
    return COMPLEMENT[first] == second


@njit(cache=True)
def _initiation_penalty_dna(base: str):
    return INITIATION_AT if (base == _A or base == _T) else INITIATION_GC


@njit(cache=True)
def _is_deoxy(sugar_index):
    return sugar_index == DEOXY


@njit(cache=True)
def hairpin(seq: np.ndarray, length: int, sugars: np.ndarray, match_dg: np.ndarray) -> tuple[float, int, int]:
    """
    Find the best step-loop structure within a single ASO.

    * A stem cannot contain a mismatch

    Returns (dG, length of stem, # of modified nucleotides within stem)
    """
    best_energy, best_stem, best_modified = 0.0, 0, 0
    for i in range(length):
        for j in range(i + MIN_LOOP + 1, length):
            if not _is_paired(seq[i], seq[j]):
                continue
            loop = j - i - 1
            if loop < MIN_LOOP or loop > MAX_LOOP:
                continue

            # The loop penalty is chemistry-blind; see santalucia.
            base_score = HAIRPIN_LOOP_PENALTY[loop]
            stack = 0.0
            modified = (0 if _is_deoxy(sugars[i]) else 1) + (0 if _is_deoxy(sugars[j]) else 1)
            for k in range(1, length):
                inner, outer = i - k, j + k
                if inner < 0 or outer >= length or not _is_paired(seq[inner], seq[outer]):
                    break
                stack += match_dg[
                    sugars[inner], sugars[inner + 1], sugars[outer], sugars[outer - 1], seq[inner], seq[inner + 1]
                ]
                modified += (0 if _is_deoxy(sugars[inner]) else 1) + (0 if _is_deoxy(sugars[outer]) else 1)

                energy = base_score + stack + _initiation_penalty_dna(seq[inner])
                if energy < best_energy:
                    best_energy, best_stem, best_modified = energy, k + 1, modified
    return best_energy, best_stem, best_modified


@njit(cache=True)
def homodimer(seq, length, sugars, match_dg):
    """
    Find the best WC ASO homodimer. Mismatches are forbidden.
    """
    best, best_run, best_modified = 0.0, 0, 0
    for shift in range(-(length - 1), length):
        start = 0
        while start < length:
            partner = length - 1 - start - shift
            if not (0 <= partner < length) or not _is_paired(seq[start], seq[partner]):
                start += 1
                continue
            run = 0.0
            end = start
            while end + 1 < length:
                next_partner = length - 2 - end - shift
                if not (0 <= next_partner < length) or not _is_paired(seq[end + 1], seq[next_partner]):
                    break
                here = length - 1 - end - shift
                run += match_dg[sugars[end], sugars[end + 1], sugars[here], sugars[here - 1], seq[end], seq[end + 1]]
                end += 1
                energy = _initiation_penalty_dna(seq[start]) + _initiation_penalty_dna(seq[end]) + run
                if energy < best:
                    modified = 0
                    for pos in range(start, end + 1):
                        against = length - 1 - pos - shift
                        modified += (0 if _is_deoxy(sugars[pos]) else 1) + (0 if _is_deoxy(sugars[against]) else 1)
                    best, best_run, best_modified = energy, end - start + 1, modified
            start += 1
    return best, best_run, best_modified


@njit(cache=True)
def _overlap_single_offset(seq, length, sugars, shift, step_dg, mismatch_dg, tandem_dg):
    """Energy of the whole overlapping region at one register, mismatches included."""

    # Setup: Check if legal overlap possible (1 nt not enough)
    first, last = -1, -1
    for i in range(length):
        if 0 <= length - 1 - i - shift < length:
            if first < 0:
                first = i
            last = i
    if first < 0 or last <= first:
        return np.nan, 0, 0

    # Now we know overlap is legal, we score it
    energy = _initiation_penalty_dna(seq[first]) + _initiation_penalty_dna(seq[last])
    crossed = 0

    if not _is_paired(seq[first], seq[length - 1 - first - shift]):
        crossed += 1

    for pos in range(first, last):
        partner1 = length - 1 - pos - shift
        partner2 = length - 2 - pos - shift
        if _is_paired(seq[pos], seq[partner1]) and _is_paired(seq[pos + 1], seq[partner2]):
            energy += step_dg[sugars[pos], sugars[pos + 1], sugars[partner1], sugars[partner2], seq[pos], seq[pos + 1]]
        elif _is_paired(seq[pos + 1], seq[partner2]):  # mismatch in the first
            energy += mismatch_dg[sugars[pos], sugars[partner1], seq[pos + 1], seq[partner2], seq[pos], seq[partner1]]
        elif _is_paired(seq[pos], seq[partner1]):  # mismatch in the second
            energy += mismatch_dg[
                sugars[pos + 1], sugars[partner2], seq[pos], seq[partner1], seq[pos + 1], seq[partner2]
            ]
            crossed += 1  # we count only mismatch in second to not double-count with mismatch in first
        else:
            # 2 mismatches don't have special weights, they get a flat penalty
            energy += tandem_dg[sugars[pos + 1], sugars[partner2]]
            crossed += 1

    return energy, last - first + 1, crossed


@njit(cache=True)
def whole_overlap(seq, length, sugars, match_dg, mismatch_dg, tandem_dg):
    """Best whole-overlap register, and the fully overlapped one reported separately.

    Returns (dG, overlap length, mismatches, fully-overlapped dG, its mismatches).
    """
    best, best_len, best_mm = np.inf, 0, 0
    flush, flush_mm = np.nan, 0
    for shift in range(-(length - 1), length):
        energy, span, crossed = _overlap_single_offset(seq, length, sugars, shift, match_dg, mismatch_dg, tandem_dg)
        if np.isnan(energy):
            continue
        if shift == 0:
            flush, flush_mm = energy, crossed
        if energy < best:
            best, best_len, best_mm = energy, span, crossed
    if not np.isfinite(best):
        return np.nan, 0, 0, flush, flush_mm
    return best, best_len, best_mm, flush, flush_mm


@njit(cache=True)
def _run_hairpin(sequences, lengths, sugars, step_dg):
    rows = sequences.shape[0]
    dg, stem, modified = np.zeros(rows), np.zeros(rows), np.zeros(rows)
    for row in range(rows):
        dg[row], stem[row], modified[row] = hairpin(sequences[row], lengths[row], sugars[row], step_dg)
    return dg, stem, modified


@njit(cache=True)
def _run_homodimer_matched(sequences, lengths, sugars, step_dg):
    rows = sequences.shape[0]
    dg, helix, modified = np.zeros(rows), np.zeros(rows), np.zeros(rows)
    for row in range(rows):
        dg[row], helix[row], modified[row] = homodimer(sequences[row], lengths[row], sugars[row], step_dg)
    return dg, helix, modified


@njit(cache=True)
def _run_homodimer_mismatched(sequences, lengths, sugars, step_dg, mismatch_dg, tandem_dg):
    rows = sequences.shape[0]
    dg, span, mismatches = np.zeros(rows), np.zeros(rows), np.zeros(rows)
    flush_dg, flush_mismatches = np.zeros(rows), np.zeros(rows)
    for row in range(rows):
        dg[row], span[row], mismatches[row], flush_dg[row], flush_mismatches[row] = whole_overlap(
            sequences[row], lengths[row], sugars[row], step_dg, mismatch_dg, tandem_dg
        )
    return dg, span, mismatches, flush_dg, flush_mismatches


def encode(sequences, patterns):
    """Sequences and their sugars as integer arrays, padded to the longest oligo."""
    index = {base: i for i, base in enumerate(DNA_BASES)}
    count = len(sequences)
    width = max(len(str(s)) for s in sequences)
    encoded = np.zeros((count, width), dtype=np.int8)
    sugars = np.zeros((count, width), dtype=np.int8)
    lengths = np.zeros(count, dtype=np.int64)
    for row, (sequence, pattern) in enumerate(zip(sequences, patterns)):
        sequence = normalize_dna(str(sequence))
        lengths[row] = len(sequence)
        for position, base in enumerate(sequence):
            encoded[row, position] = index[base]
        pattern = str(pattern)
        for position in range(min(len(pattern), len(sequence))):
            sugars[row, position] = SUGARS.get(pattern[position].upper(), DEOXY)
    return encoded, lengths, sugars


def calculate_self_aso(sequences, patterns):
    """Every self-structure quantity for each oligo, keyed by name."""
    encoded, lengths, sugars = encode(sequences, patterns)
    hp_dg, hp_stem, hp_modpairs = _run_hairpin(encoded, lengths, sugars, DINUCLEOTIDE_ENERGY)
    homodimer_perfect_dg, homodimer_perfect_helix, homodimer_perfect_modpairs = _run_homodimer_matched(
        encoded, lengths, sugars, DINUCLEOTIDE_ENERGY
    )
    homodimer_mismatched_dg, homodimer_mismatched_len, homodimer_mismatched_nmm, full_dg, full_nmm = (
        _run_homodimer_mismatched(
            encoded, lengths, sugars, DINUCLEOTIDE_ENERGY, MISMATCH_ENERGY, TANDEM_MISMATCH_ENERGY
        )
    )
    return {
        "hp_dg": hp_dg,
        "hp_stem": hp_stem,
        "hp_modpairs": hp_modpairs,
        "hp_dg_per_nt": hp_dg / lengths,
        "homodimer_perfect_dg": homodimer_perfect_dg,
        "homodimer_perfect_helix": homodimer_perfect_helix,
        "homodimer_perfect_modpairs": homodimer_perfect_modpairs,
        "homodimer_perfect_dg_per_nt": homodimer_perfect_dg / lengths,
        "homodimer_perfect_minus_hp": homodimer_perfect_dg - hp_dg,
        "homodimer_mismatched_dg": homodimer_mismatched_dg,
        "homodimer_mismatched_len": homodimer_mismatched_len,
        "homodimer_mismatched_nmm": homodimer_mismatched_nmm,
        "homodimer_mismatched_dg_per_nt": homodimer_mismatched_dg / lengths,
        "homodimer_mismatched_full_dg": full_dg,
        "homodimer_mismatched_full_nmm": full_nmm,
    }
