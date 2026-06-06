"""2'-MOE nearest-neighbour weights from MD simulations of 2'-MOE-modified DNA/RNA hybrid
duplexes (the MOE sits on the DNA antisense strand, paired against the RNA target -- the
gapmer wing-vs-target context).

Source: Park et al., "Thermodynamic Parameter Estimation for Modified Oligonucleotides Using
Molecular Dynamics Simulations," J. Phys. Chem. B 2025, 129(11):2934-2945
(doi:10.1021/acs.jpcb.4c08344). The table lives next to this module as weights/moe_md.csv.
Each dinucleotide key carries a sugar marker per base: ' = unmodified, * = 2'-MOE. Mixed-prime
keys cover the MOE/unmodified junction at the boundary between a MOE wing and the DNA gap.
"""

from importlib.resources import files

import pandas as pd

from ...common.modifications import get_longest_dna_gap
from ...util import celsius_to_kelvin

_BODY_TEMPERATURE_K = celsius_to_kelvin(37.0)
_WEIGHTS = files("tauso.features.hybridization") / "weights"


def _gibbs_37_weights(table, h_col, s_col):
    weights = {}
    for _, row in table.iterrows():
        weights[row["nucleotide"]] = row[h_col] - _BODY_TEMPERATURE_K * (row[s_col] / 1000.0)
    return weights


with (_WEIGHTS / "moe_md.csv").open() as _handle:
    _moe_md = pd.read_csv(_handle, comment="#")

MOE_MD_PB_WEIGHTS = _gibbs_37_weights(_moe_md, "H_PB", "S_PB")
MOE_MD_GB_WEIGHTS = _gibbs_37_weights(_moe_md, "H_GB", "S_GB")


def _moe_weights(simul_type: str) -> dict:
    if simul_type == "gb":
        return MOE_MD_GB_WEIGHTS
    if simul_type == "pb":
        return MOE_MD_PB_WEIGHTS
    raise ValueError(f"Unknown simulation type: {simul_type}")


def get_moe_md_contribution(seq: str, chemical_pattern, modification, simul_type="gb", region=None):
    """Whole-oligo 2'-MOE-modified DNA/RNA affinity from MD (kcal/mol): the MD nearest-neighbour
    energy summed over the full duplex (MOE-bearing stacks plus the unmodified DNA-gap stacks),
    so it is the modified ASO:RNA affinity rather than a wing-only increment. Note this is a
    fully MD-derived quantity (gap included), not on the same footing as the physical DNA/RNA
    weights. NaN for non-MOE oligos. ``region`` restricts the sum to one wing ('wing5'/'wing3').
    """
    if "MOE" not in modification:
        return float("nan")
    if not isinstance(chemical_pattern, str) or len(chemical_pattern) != len(seq):
        return float("nan")

    seq = seq.replace("T", "U")
    weights = _moe_weights(simul_type)
    gap_start, gap_end, gap_len = get_longest_dna_gap(chemical_pattern)
    if gap_len == 0:
        gap_start, gap_end = len(seq), len(seq)

    total = 0.0
    for i in range(len(seq) - 1):
        if region == "wing5" and not i < gap_start:
            continue
        if region == "wing3" and not i >= gap_end:
            continue
        L = seq[i] + ("*" if chemical_pattern[i] == "M" else "'")
        R = seq[i + 1] + ("*" if chemical_pattern[i + 1] == "M" else "'")
        total += weights[L + R]
    return total
