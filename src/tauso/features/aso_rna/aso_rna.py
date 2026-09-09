"""Duplex geometry of the ASO against its RNA target, averaged over each region of the gapmer.

A gapmer wing is 2'-MOE against RNA and its gap is DNA against RNA, so every step of the
duplex is looked up in a table built for the sugar pair it actually has. Steps where the
chemistry changes use the junction cells, which are keyed on both sugars.

Thirteen geometric observables are read at every step; each is averaged over the 5' wing, the
gap and the 3' wing. The averaging is done twice, once over the tables' mean values and once
over their spread, since how variable a step's geometry is carries information the mean does
not.

Steps straddling a region boundary belong to no region, so a value is averaged only over
steps that sit wholly inside one.
"""

import csv
from importlib import resources

import numpy as np

from ...common.modifications import get_longest_dna_gap
from ...util import normalize_dna

OBSERVABLES = (
    "Roll",
    "Tilt",
    "hIncl",
    "hTip",
    "hY",
    "hRise",
    "Shift",
    "Rise",
    "Slide",
    "Zp",
    "Twist",
    "hTwist",
    "hX",
)
REGIONS = ("wing5", "gap", "wing3")

SUGARS = {"M": "M", "C": "E", "O": "O", "R": "R"}
"""`chemical_pattern` letters mapped to the sugar the weight tables are keyed on.

The column writes cEt as "C" where the tables name it "E". Anything absent is deoxy.
"""

FEATURE_NAMES = [f"{o}_{r}" for o in OBSERVABLES for r in REGIONS] + [
    f"{o}_{r}_spread" for o in OBSERVABLES for r in REGIONS
]


def _load(filename, stat):
    tables = {o: {} for o in OBSERVABLES}
    with resources.files(__package__).joinpath("weights", filename).open() as handle:
        for row in csv.DictReader(handle):
            if row["stat"] == stat and row["obs"] in tables:
                tables[row["obs"]][row["cell"]] = float(row["value"])
    return tables


UNIFORM_MEAN = _load("rna_uniform.csv", "mean")
UNIFORM_SPREAD = _load("rna_uniform.csv", "sd")
JUNCTION_MEAN = _load("rna_junction.csv", "mean")


def sugars(chemical_pattern, length):
    text = str(chemical_pattern)
    return [SUGARS.get(text[i].upper(), "D") if i < len(text) else "D" for i in range(length)]


def step_regions(chemical_pattern, length):
    """Boolean mask per region over the L-1 steps; a boundary step is in none of them."""
    start, end, found = get_longest_dna_gap(str(chemical_pattern))
    if not found:
        start, end = 0, length
    steps = np.arange(1, length)
    return {
        "wing5": (steps - 1 < start) & (steps < start),
        "gap": (steps - 1 >= start) & (steps < end),
        "wing3": (steps - 1 >= end) & (steps >= end),
    }


def step_cell(sugar_5, sugar_3, dinucleotide):
    """Which table cell a step is scored in, given the sugars on either side of it."""
    if sugar_5 == sugar_3:
        return f"{sugar_5}:R@{dinucleotide}", False
    return f"{sugar_5}{sugar_3}|RR@{dinucleotide}", True


def _profile(sequence, sugar, uniform):
    """One value per step for every observable, or NaN where the cell is unmeasured."""
    length = len(sequence)
    values = {o: np.full(length - 1, np.nan) for o in OBSERVABLES}
    for i in range(1, length):
        cell, is_junction = step_cell(sugar[i - 1], sugar[i], sequence[i - 1] + sequence[i])
        tables = JUNCTION_MEAN if is_junction else uniform
        for o in OBSERVABLES:
            values[o][i - 1] = tables[o].get(cell, np.nan)
    return values


def calculate_aso_rna(sequences, chemical_patterns):
    """One array per feature name, each as long as `sequences`."""
    if len(sequences) != len(chemical_patterns):
        raise ValueError(f"got {len(sequences)} sequences against {len(chemical_patterns)} chemical patterns")

    out = {name: np.full(len(sequences), np.nan) for name in FEATURE_NAMES}
    for row, (sequence, pattern) in enumerate(zip(sequences, chemical_patterns)):
        sequence = normalize_dna(str(sequence))
        length = len(sequence)
        if length < 2:
            continue
        sugar = sugars(pattern, length)
        masks = step_regions(pattern, length)
        for uniform, suffix in ((UNIFORM_MEAN, ""), (UNIFORM_SPREAD, "_spread")):
            values = _profile(sequence, sugar, uniform)
            for o in OBSERVABLES:
                for region, mask in masks.items():
                    inside = values[o][mask]
                    inside = inside[np.isfinite(inside)]
                    if inside.size:
                        out[f"{o}_{region}{suffix}"][row] = inside.mean()
    return out
