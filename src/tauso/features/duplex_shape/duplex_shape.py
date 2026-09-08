"""Duplex geometry along the oligo, predicted from sequence by nearest-neighbour lookup.

`weights/duplex_shape.csv` gives the mean of an observable for each dinucleotide in each
chemistry state, measured on 480 phosphodiester duplexes. Applying it to a sequence gives a
value per step; reducing that to a few numbers per oligo is what the features are.

Thirteen observables are carried, in fourteen columns -- Slide contributes both its mean and
its fluctuation. They are the ones that predict held-out duplexes in every chemistry state.

The reduction keeps the two wings apart rather than averaging them. A gapmer's ends are not
interchangeable, and pooling them was measured to lose most of what the geometry contributes.
A step belongs to a region only if both of the positions it spans do, so the two steps that
straddle a wing/gap boundary count towards neither.
"""

from importlib.resources import files

import numpy as np
import pandas as pd

from ...common.modifications import get_longest_dna_gap
from ...util import normalize_dna

_WEIGHTS = files("tauso.features.duplex_shape") / "weights"

SUGAR_STATE = {"M": "M", "C": "E"}
"""`chemical_pattern` letters, as the sugar codes the weight table is keyed on. Deoxy is
anything else, which is how the column writes it."""

REGIONS = ("wing5", "gap", "wing3")

OBSERVABLES = (
    "roll",
    "hincl",
    "htip",
    "tilt",
    "hrise",
    "hy",
    "shift",
    "slide",
    "propeller",
    "rise",
    "buckle",
    "stagger",
    "zp",
    "slidesd",
)

_COLUMNS = {
    "roll": ("Roll", "mean"),
    "hincl": ("hIncl", "mean"),
    "htip": ("hTip", "mean"),
    "tilt": ("Tilt", "mean"),
    "hrise": ("hRise", "mean"),
    "hy": ("hY", "mean"),
    "shift": ("Shift", "mean"),
    "slide": ("Slide", "mean"),
    "propeller": ("Propeller", "mean"),
    "rise": ("Rise", "mean"),
    "buckle": ("Buckle", "mean"),
    "stagger": ("Stagger", "mean"),
    "zp": ("Zp", "mean"),
    "slidesd": ("Slide", "sd"),
}

QUANTITIES = [f"{observable}_{region}" for observable in OBSERVABLES for region in REGIONS]


def _read():
    with (_WEIGHTS / "duplex_shape.csv").open() as handle:
        return pd.read_csv(handle)


def build_tables():
    """One table per observable, keyed by (chemistry state, dinucleotide)."""
    table = _read()
    tables = {}
    for name, (observable, statistic) in _COLUMNS.items():
        rows = table[(table.obs == observable) & (table.stat == statistic)]
        tables[name] = {(row.state, row.cell): row.value for row in rows.itertuples()}
    return tables


TABLES = build_tables()


def sugars(pattern, length):
    text = str(pattern)
    return [SUGAR_STATE.get(text[i].upper(), "D") if i < len(text) else "D" for i in range(length)]


def step_masks(pattern, length):
    """Which steps lie inside each region. A step spans two positions and needs both."""
    start, end, found = get_longest_dna_gap(str(pattern))
    if not found:
        start, end = 0, length
    steps = np.arange(1, length)
    return {
        "wing5": (steps - 1 < start) & (steps < start),
        "gap": (steps - 1 >= start) & (steps < end),
        "wing3": (steps - 1 >= end) & (steps >= end),
    }


def score_one(sequence, pattern):
    """Every region mean for one oligo, keyed by quantity. NaN where a region has no step."""
    sequence = normalize_dna(str(sequence))
    length = len(sequence)
    sugar = sugars(pattern, length)
    masks = step_masks(pattern, length)
    dinucleotides = [sequence[i - 1] + sequence[i] for i in range(1, length)]
    states = [sorted({f"{sugar[i - 1]}:{sugar[i - 1]}", f"{sugar[i]}:{sugar[i]}"}) for i in range(1, length)]

    scored = {}
    for name in OBSERVABLES:
        table = TABLES[name]
        values = np.array(
            [
                np.mean([v for v in (table.get((s, d)) for s in state) if v is not None] or [np.nan])
                for d, state in zip(dinucleotides, states)
            ]
        )
        for region, mask in masks.items():
            selected = values[mask]
            scored[f"{name}_{region}"] = float(np.nanmean(selected)) if selected.size else np.nan
    return scored


def calculate_duplex_shape(sequences, patterns):
    """Every duplex-shape quantity for each oligo, keyed by name."""
    scored = [score_one(sequence, pattern) for sequence, pattern in zip(sequences, patterns)]
    return {quantity: np.array([row[quantity] for row in scored]) for quantity in QUANTITIES}
