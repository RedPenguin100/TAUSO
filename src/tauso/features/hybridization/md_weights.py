"""2'-MOE nearest-neighbour weights from MD simulations of MOE-modified RNA duplexes.

Source: https://pubs.acs.org/doi/10.1021/acs.jpcb.4c08344
The table lives next to this module as weights/moe_md.csv. Each dinucleotide key
carries a sugar marker per base: ' = unmodified, * = 2'-MOE. Mixed-prime keys cover
the MOE/unmodified junction at the boundary between a MOE wing and the DNA gap.
"""

from importlib.resources import files

import pandas as pd

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


def get_moe_md_contribution(seq: str, chemical_pattern, modification, simul_type="gb"):
    """Sum of MD nearest-neighbour weights over the 2'-MOE-bearing dinucleotides (kcal/mol).

    Only defined for MOE oligos; returns 0 otherwise. Dinucleotides where both sugars
    are unmodified are skipped, so this isolates the contribution of the 2'-MOE wings
    and their junctions with the central DNA gap.
    """
    if "MOE" not in modification:
        return 0

    seq = seq.replace("T", "U")
    weights = _moe_weights(simul_type)

    total = 0.0
    for i in range(len(seq) - 1):
        L = seq[i] + ("*" if chemical_pattern[i] == "M" else "'")
        R = seq[i + 1] + ("*" if chemical_pattern[i + 1] == "M" else "'")
        if "'" in L and "'" in R:
            continue
        total += weights[L + R]
    return total
