"""Nearest-neighbour weights over three sugar chemistries: deoxy, 2'-MOE and cEt.

`nn_table_all.csv` carries the same-chemistry and cross classes (DD, MM, MD, EE, ED) and
`nn_junction_by_step.csv` the fourteen junction subtypes, which stand in for a single lumped
JUNCTION column. Both are fitted on phosphorothioate steps and calibrated on PO-DNA against
SantaLucia, so they already describe a PS duplex on the published kcal/mol scale.

A step is classified from the sugars of all four nucleotides it spans:

    all four alike             DD / MM / EE
    each strand alike, differ  MD / ED, the dinucleotide flipped when the modified strand is
                               the second one, matching how the fit pooled the orientations
    otherwise                  the junction subtype for that pattern

Junction patterns canonicalise under ``a1a2|b1b2 -> b2b1|a2a1`` with the dinucleotide taken to
its reverse complement: the same duplex read from the other strand.

Initiation and the hairpin loop penalty stay at SantaLucia, since the fit contains neither
helix ends free of fraying nor any loops.
"""

from importlib.resources import files

import numpy as np
import pandas as pd

from ...util import DNA_BASES, WATSON_CRICK_MAP

DEOXY, MOE, CET = 0, 1, 2
SUGAR_CODE = {DEOXY: "D", MOE: "M", CET: "E"}

COMPLEMENT = np.array([DNA_BASES.index(WATSON_CRICK_MAP[b]) for b in DNA_BASES], dtype=np.int8)
"""Watson-Crick partner of each base, as indices into DNA_BASES.

Taken from the pairing map in `util` so the two cannot drift apart, and kept as an array
because the kernels reading it are numba-compiled and cannot look inside a dict.
"""

INIT_GC, INIT_AT = 0.98, 1.03
"""SantaLucia helix initiation, kcal/mol."""

HAIRPIN_LOOP = np.array(
    [
        0.0,
        0.0,
        0.0,
        3.5,
        3.5,
        3.3,
        4.0,
        4.2,
        4.3,
        4.5,
        4.6,
        4.7,
        4.8,
        4.9,
        4.9,
        5.0,
        5.1,
        5.1,
        5.2,
        5.2,
        5.3,
        5.3,
        5.4,
        5.4,
        5.4,
        5.5,
        5.5,
        5.5,
        5.6,
        5.6,
        5.6,
    ],
    dtype=np.float64,
)
"""Loop penalty by loop length, kcal/mol. Index 0-2 are unreachable: a loop needs three bases."""

MIN_LOOP, MAX_LOOP = 3, len(HAIRPIN_LOOP) - 1

_WEIGHTS = files("tauso.features.self_structure") / "weights"
_BASE_CLASSES = ["DD", "MM", "MD", "EE", "ED", "JUNCTION"]


def _read(name):
    with (_WEIGHTS / name).open() as handle:
        return pd.read_csv(handle)


def build_tables():
    """Return (parameters[class, base, base], class_map[a1, a2, b1, b2] -> (class, flip)).

    `flip` marks the steps whose dinucleotide must be read as its reverse complement, which is
    how the fit pooled the two strand orientations of an asymmetric class.
    """
    table = _read("nn_table_all.csv").set_index("din")
    junctions = _read("nn_junction_by_step.csv")

    # The fit produced no GC measurement for the cross classes; the mean of the two
    # same-chemistry classes stands in.
    for cross, same in (("MD", "MM"), ("ED", "EE")):
        missing = table[cross].isna()
        if missing.any():
            table.loc[missing, cross] = ((table["DD"] + table[same]) / 2)[missing]

    patterns = sorted(set(zip(junctions.chem, junctions.pattern)))
    names = _BASE_CLASSES + [f"{chem}:{pattern}" for chem, pattern in patterns]
    parameters = np.zeros((len(names), 4, 4), dtype=np.float64)
    for k, cls in enumerate(_BASE_CLASSES):
        for din in table.index:
            parameters[k, DNA_BASES.index(din[0]), DNA_BASES.index(din[1])] = table[cls][din]
    lumped = table["JUNCTION"]
    for k, (chem, pattern) in enumerate(patterns, start=len(_BASE_CLASSES)):
        sub = junctions[(junctions.chem == chem) & (junctions.pattern == pattern)].set_index("din")
        for din in table.index:
            value = sub.dG[din] if din in sub.index else lumped[din]
            parameters[k, DNA_BASES.index(din[0]), DNA_BASES.index(din[1])] = value

    index = {name: i for i, name in enumerate(names)}
    class_map = np.zeros((3, 3, 3, 3, 2), dtype=np.int64)
    for a1 in (DEOXY, MOE, CET):
        for a2 in (DEOXY, MOE, CET):
            for b1 in (DEOXY, MOE, CET):
                for b2 in (DEOXY, MOE, CET):
                    class_map[a1, a2, b1, b2] = _classify(a1, a2, b1, b2, index)
    return parameters, class_map


def _classify(a1, a2, b1, b2, index):
    if a1 == a2 == b1 == b2:
        return index[{DEOXY: "DD", MOE: "MM", CET: "EE"}[a1]], 0
    if a1 == a2 and b1 == b2:
        modified = a1 if a1 != DEOXY else b1
        if modified == DEOXY:
            return index["DD"], 0
        return index["MD" if modified == MOE else "ED"], 0 if a1 != DEOXY else 1
    sugars = tuple(SUGAR_CODE[x] for x in (a1, a2, b1, b2))
    chem = "M" if MOE in (a1, a2, b1, b2) else "E"
    forward = "{}{}|{}{}".format(*sugars)
    reverse = f"{sugars[3]}{sugars[2]}|{sugars[1]}{sugars[0]}"
    if f"{chem}:{forward}" in index:
        return index[f"{chem}:{forward}"], 0
    if f"{chem}:{reverse}" in index:
        return index[f"{chem}:{reverse}"], 1
    # Mixed 2'-MOE and cEt in one step, which the fit does not cover.
    return index["JUNCTION"], 0


PARAMETERS, CLASS_MAP = build_tables()
