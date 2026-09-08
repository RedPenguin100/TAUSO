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

from ...common.modifications import check_pattern_length, get_longest_dna_gap
from ...util import celsius_to_kelvin, dna_to_rna

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


# The weight keys interleave base and sugar marker per residue, while the sum walks a sequence
# and a chemical pattern side by side. Re-keying by the two slices it already holds -- the RNA
# dinucleotide followed by its two pattern characters -- makes a stack one dict lookup.
_MOE_BY_DINUCLEOTIDE = {
    simul: {
        b1 + b2 + c1 + c2: weights[b1 + ("*" if c1 == "M" else "'") + b2 + ("*" if c2 == "M" else "'")]
        for b1 in "ACGU"
        for b2 in "ACGU"
        for c1 in "Md"
        for c2 in "Md"
    }
    for simul, weights in (("gb", MOE_MD_GB_WEIGHTS), ("pb", MOE_MD_PB_WEIGHTS))
}


def get_moe_md_contribution(seq: str, chemical_pattern, modification, simul_type="gb", region=None):
    """Whole-oligo 2'-MOE-modified DNA/RNA affinity from MD (kcal/mol): the MD nearest-neighbour
    energy summed over the full duplex (MOE-bearing stacks plus the unmodified DNA-gap stacks),
    so it is the modified ASO:RNA affinity rather than a wing-only increment. Note this is a
    fully MD-derived quantity (gap included), not on the same footing as the physical DNA/RNA
    weights. NaN for non-MOE oligos. ``region`` restricts the sum to one wing ('wing5'/'wing3').
    """
    if "MOE" not in modification:
        return float("nan")
    check_pattern_length(seq, chemical_pattern)
    if not isinstance(chemical_pattern, str):
        return float("nan")

    seq = dna_to_rna(seq)
    weights = _moe_weights(simul_type)
    table = _MOE_BY_DINUCLEOTIDE[simul_type]
    gap_start, gap_end, gap_len = get_longest_dna_gap(chemical_pattern)
    if gap_len == 0:
        gap_start, gap_end = len(seq), len(seq)

    if region == "wing5":
        start, stop = 0, gap_start
    elif region == "wing3":
        start, stop = gap_end, len(seq)
    else:
        start, stop = 0, len(seq)

    total = 0.0
    for i in range(start, min(stop, len(seq) - 1)):
        key = seq[i : i + 2] + chemical_pattern[i : i + 2]
        if key in table:
            total += table[key]
        else:
            L = seq[i] + ("*" if chemical_pattern[i] == "M" else "'")
            R = seq[i + 1] + ("*" if chemical_pattern[i + 1] == "M" else "'")
            total += weights[L + R]
    return total
