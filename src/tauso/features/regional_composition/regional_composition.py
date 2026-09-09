"""Dinucleotide composition of each region of a gapmer, against the rest of the oligo.

The oligo is cut into five regions by the DNA gap: the 5' wing, the step that crosses into
the gap, the gap, the step that crosses out, and the 3' wing. Each region gets the frequency
of all sixteen dinucleotides among its steps, minus that dinucleotide's frequency over the
steps outside the region.

Everything is anchored to the gap boundary rather than to the oligo ends, which is what the
one-hot terminal features cannot express: wing length varies across the corpus, so a fixed
offset from the 5' end lands in different regions for different oligos.

A step spans two adjacent bases, so an oligo of length L has L-1 of them and the five regions
partition those exactly. A region is compared against the steps outside it rather than against
every step, so it does not sit inside its own reference.

NaN, rather than zero, stands for a comparison that cannot be made: a region with no steps,
such as both wings of an all-DNA oligo, and a region holding every step, which leaves nothing
to compare against.
"""

import numpy as np

from ...common.modifications import get_longest_dna_gap
from ...util import normalize_dna

BASES = "ACGT"
DINUCLEOTIDES = [a + b for a in BASES for b in BASES]
REGIONS = ("wing5", "j5", "gap", "j3", "wing3")

FEATURE_NAMES = [f"{d}_{r}" for r in REGIONS for d in DINUCLEOTIDES]


def step_regions(chemical_pattern, length):
    """Boolean mask per region over the L-1 steps, each step in exactly one region."""
    start, end, found = get_longest_dna_gap(str(chemical_pattern))
    if not found:
        start, end = 0, length
    steps = np.arange(1, length)
    return {
        "wing5": (steps - 1 < start) & (steps < start),
        "j5": steps == start,
        "gap": (steps - 1 >= start) & (steps < end),
        "j3": steps == end,
        "wing3": (steps - 1 >= end) & (steps >= end),
    }


def _counts(sequence, steps):
    """How many times each dinucleotide appears among the given steps, and how many there are."""
    counts = dict.fromkeys(DINUCLEOTIDES, 0)
    total = 0
    for i in steps:
        pair = sequence[i - 1] + sequence[i]
        if pair in counts:
            counts[pair] += 1
            total += 1
    return counts, total


def calculate_regional_composition(sequences, chemical_patterns):
    """One array per feature name, each as long as `sequences`."""
    if len(sequences) != len(chemical_patterns):
        raise ValueError(f"got {len(sequences)} sequences against {len(chemical_patterns)} chemical patterns")

    out = {name: np.full(len(sequences), np.nan) for name in FEATURE_NAMES}
    for row, (sequence, pattern) in enumerate(zip(sequences, chemical_patterns)):
        sequence = normalize_dna(str(sequence))
        length = len(sequence)
        if length < 2:
            continue
        masks = step_regions(pattern, length)
        whole, whole_total = _counts(sequence, np.arange(1, length))
        for region, mask in masks.items():
            local, local_total = _counts(sequence, np.where(mask)[0] + 1)
            outside_total = whole_total - local_total
            if not local_total or not outside_total:
                continue
            for d in DINUCLEOTIDES:
                outside = (whole[d] - local[d]) / outside_total
                out[f"{d}_{region}"][row] = local[d] / local_total - outside
    return out
