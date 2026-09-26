"""Choosing sequences to measure with.

Shared by both design steps. A sequence is fit to measure with if it is not dominated by one base,
not prone to folding on itself, and carries enough distinct dinucleotides to be worth a duplex --
`usable` below. A *set* of sequences is fit if it covers the sixteen dinucleotides evenly, which
is what `balanced` arranges.

`usable` describes a sequence as it will be simulated, so any step that writes into a sequence
must run it again afterwards: a filter passed before the write says nothing about what comes out.
"""
from collections import Counter

import numpy as np

from consts import DUPLEX_LENGTH as L, reverse_complement

MARGIN = 2                  # steps this close to either end fray and are not measured
BASES = list("ACGT")
DINUCLEOTIDES = [a + b for a in "ACGT" for b in "ACGT"]
SEQUENCE_SEED = 50          # the campaign this design belongs to
N_SEQUENCES = 48
CANDIDATE_POOL = 600        # drawn, then thinned to the set that covers the cells most evenly


def self_complementarity(sequence):
    """Longest stretch that would pair with the strand's own reverse complement."""
    other = reverse_complement(sequence)
    best = 0
    for offset in range(-(L - 4), L - 3):
        run = 0
        for i in range(L):
            j = i + offset
            if 0 <= j < L and sequence[i] == other[j]:
                run += 1
                best = max(best, run)
            else:
                run = 0
    return best


def has_tetranucleotide(sequence):
    """Whether any base is repeated four times in a row."""
    return any(base * 4 in sequence for base in BASES)


def gc_content(sequence):
    return sum(b in "GC" for b in sequence) / len(sequence)


def dinucleotide_count(sequence):
    return len({sequence[i:i + 2] for i in range(L - 1)})


def usable(sequence):
    """Whether a candidate is fit to measure with: no long homopolymer, balanced GC, not prone to
    folding on itself, and carrying enough distinct dinucleotides to be worth a duplex."""
    if has_tetranucleotide(sequence):
        return False
    if not 0.31 < gc_content(sequence) < 0.69:
        return False
    if self_complementarity(sequence) >= 6:
        return False
    return dinucleotide_count(sequence) >= 10


def make_sequences(n, rng):
    """`n` distinct sequences passing `usable`, drawn from `rng`."""
    out = []
    while len(out) < n:
        candidate = "".join(rng.choice(BASES, L))
        if usable(candidate) and candidate not in out:
            out.append(candidate)
    return out


def interior_steps(sequence):
    """The dinucleotides a duplex actually contributes, terminal steps excluded."""
    return [sequence[p - 1:p + 1] for p in range(MARGIN + 1, L - MARGIN)]


def unevenness(coverage):
    """How lopsided a coverage count is. Lower is better; ties go to the fullest thin cell."""
    counts = np.array([coverage[c] for c in DINUCLEOTIDES])
    return counts.std(), -counts.min()


def balanced(pool, n):
    """Pick `n` sequences from `pool` so the sixteen dinucleotides are measured equally often.

    A random sequence does not know which cells are short, so taking the first `n` that pass
    `usable` leaves some cells at 21 steps and others at 42. Picking greedily -- each round, keep
    whichever candidate leaves coverage least lopsided -- brings every cell to 32-34 for the same
    number of duplexes, and every sequence is still an ordinary one drawn under the filter.

    Whole sequences are compared rather than single bases because steps overlap: the base at one
    position closes one step and opens the next, so a base-by-base rule cannot see what it is
    about to cost.
    """
    chosen = []
    coverage = Counter()
    remaining = list(pool)

    for _ in range(n):
        def coverage_with(candidate):
            return unevenness(coverage + Counter(interior_steps(candidate)))

        pick = min(remaining, key=coverage_with)
        remaining.remove(pick)
        chosen.append(pick)
        coverage.update(interior_steps(pick))

    return chosen
