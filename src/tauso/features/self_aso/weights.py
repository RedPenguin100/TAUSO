"""Nearest-neighbour weights for the self-structure search, resolved into flat lookup tables.

Stacking is SantaLucia's unified DNA/DNA set at 37 C. Each modified sugar then adds a constant
bonus, applied per modified nucleotide the step spans. A nucleotide sits in two steps, so the
coefficient carried per slot is half the per-modification value.

The bonuses are on the scale of the published melting measurements -- 2'-MOE from Freier &
Altmann, cEt inferred from LNA -- and were chosen from a sweep of the two over model accuracy.
They stand in for sugar-specific nearest-neighbour tables, which do not exist for cEt at all.

Everything here is resolved once at import so the hot loop indexes arrays directly.
"""

import numpy as np

from ...util import DNA_BASES, WATSON_CRICK_MAP
from ..hybridization.weights.dna import DNA_DNA_WEIGHTS

DEOXY, MOE, CET = 0, 1, 2

BODY_TEMPERATURE_K = 310.15

MODIFICATION_BONUS = {MOE: 0.8, CET: 2.5}
"""Stabilisation each modified nucleotide adds, kcal/mol."""

MISMATCH_PENALTY = 0.4
"""Cost of carrying one mismatched pair through a helix, kcal/mol.

A working constant, not a measured value: no mismatch set covers these sugars. It replaces a
per-dinucleotide table that was measured to add nothing over a constant.
"""

COMPLEMENT = np.array([DNA_BASES.index(WATSON_CRICK_MAP[b]) for b in DNA_BASES], dtype=np.int8)
"""Watson-Crick partner of each base, as indices into DNA_BASES.

Taken from the pairing map in `util` so the two cannot drift apart, and kept as an array
because the kernels reading it are numba-compiled and cannot look inside a dict.
"""


def _slot_bonus(sugar):
    """Half the per-modification bonus, since each nucleotide is spanned by two steps."""
    return -MODIFICATION_BONUS.get(sugar, 0.0) / 2.0


def build_dna_plane():
    """SantaLucia dG37 for every dinucleotide, kcal/mol, indexed [base, base]."""
    plane = np.zeros((4, 4), dtype=np.float64)
    for first in DNA_BASES:
        for second in DNA_BASES:
            key = f"{first}{second}/{WATSON_CRICK_MAP[first]}{WATSON_CRICK_MAP[second]}"
            entry = DNA_DNA_WEIGHTS[key]
            plane[DNA_BASES.index(first), DNA_BASES.index(second)] = (
                entry["dH"] - BODY_TEMPERATURE_K * entry["dS"] / 1000.0
            )
    return plane


DNA_PLANE = build_dna_plane()


def build_step_energy():
    """Stacking energy for every sugar and base combination, kcal/mol.

    Indexed `[top1, top2, bottom1, bottom2, base, base]` over the four sugars the step spans.
    """
    energy = np.zeros((3, 3, 3, 3, 4, 4), dtype=np.float64)
    for top1 in range(3):
        for top2 in range(3):
            for bottom1 in range(3):
                for bottom2 in range(3):
                    bonus = sum(_slot_bonus(s) for s in (top1, top2, bottom1, bottom2))
                    energy[top1, top2, bottom1, bottom2] = DNA_PLANE + bonus
    return energy


def build_mismatch_energy():
    """Cost of crossing one mismatched pair, kcal/mol.

    Indexed `[top sugar, bottom sugar, anchor base, anchor base, mismatch base, mismatch base]`.
    The modification bonus is a Watson-Crick stacking term and is not credited here: a
    mismatched pair is not stacked, and crediting it would let the search cross mismatches to
    collect the bonus.
    """
    return np.full((3, 3, 4, 4, 4, 4), MISMATCH_PENALTY, dtype=np.float64)


def build_tandem_energy():
    """Cost of a mismatch with no Watson-Crick pair beside it, kcal/mol, `[top, bottom]`.

    Two mismatches in a row are not two independent stacks, so this takes no bases, and like
    the anchored table it carries no modification bonus.
    """
    return np.full((3, 3), MISMATCH_PENALTY, dtype=np.float64)


DINUCLEOTIDE_ENERGY = build_step_energy()
MISMATCH_ENERGY = build_mismatch_energy()
TANDEM_MISMATCH_ENERGY = build_tandem_energy()
