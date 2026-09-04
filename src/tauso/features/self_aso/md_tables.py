"""The TAUSO MD nearest-neighbour fit, resolved into flat lookup tables.

The CSVs in `tauso_md_weights` are the fit's own output. Everything here reads them once at
import and resolves them into arrays the search indexes directly, so no decision about which
fit covers a step, or which orientation to read it in, is left to the hot loop.

Two fits are carried. `nn_weights_final.csv` is the later one and is preferred; it names the
sugar and linkage of both strands of a step, so a deoxy strand paired against 2'-MOE has its
own cell. The older pair of tables covers what the later fit does not -- cEt against deoxy or
2'-MOE, and any step whose two positions on one strand differ in sugar -- and is described
below.

Steps are scored as phosphodiester throughout. The later fit does separate the linkages, but
the difference is small (~0.1 kcal/mol) and measured to change nothing downstream, while
reading everything from the PO cells reaches the fitted table on 70% of steps rather than 62%:
the fit has more cross-chemistry PO cells than it has backbone-specific ones.

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

Everything the searches need that this fit does not measure -- helix initiation and the
hairpin loop penalty -- lives in `santalucia`.
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

_WEIGHTS = files("tauso.features.self_aso") / "tauso_md_weights"
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


# --- the later fit -------------------------------------------------------------------

_WC_PAIRS = ["AT", "CG", "GC", "TA"]
_FINAL_SUGAR = {"D": DEOXY, "M": MOE, "E": CET}


def _canonical_step(step):
    """A step and its reverse complement are one stack read from either strand."""
    reverse = WATSON_CRICK_MAP[step[1]] + WATSON_CRICK_MAP[step[0]]
    return min(step, reverse)


def build_final_tables():
    """Return (watson_crick[top, bottom, base, base], mismatch[top, bottom, anchor, base, base, side]).

    Both are indexed by the sugar of each strand and hold NaN where the fit has no cell, which
    is how the callers know to fall back. `anchor` is the Watson-Crick pair the mismatch leans
    on, indexed by its top base; `side` is 0 when the mismatch sits 3' of that anchor and 1
    when it sits 5'.
    """
    rows = _read("nn_weights_final.csv")
    watson_crick = np.full((3, 3, 4, 4), np.nan)
    mismatch = np.full((3, 3, 4, 4, 4, 2), np.nan)

    for row in rows.itertuples():
        left, right = (part.strip() for part in row.state.split(":"))
        (top_sugar, top_link), (bottom_sugar, bottom_link) = left.split("-"), right.split("-")
        if top_link != "PO" or bottom_link != "PO":
            continue
        top, bottom = _FINAL_SUGAR[top_sugar], _FINAL_SUGAR[bottom_sugar]
        for a, b in ((top, bottom), (bottom, top)):
            if row.kind == "wc":
                for first in DNA_BASES:
                    for second in DNA_BASES:
                        if _canonical_step(first + second) == row.label:
                            watson_crick[a, b, DNA_BASES.index(first), DNA_BASES.index(second)] = row.dG
            else:
                anchor, pair, side = row.label.split("|")
                mismatch[
                    a,
                    b,
                    _WC_PAIRS.index(anchor),
                    DNA_BASES.index(pair[0]),
                    DNA_BASES.index(pair[1]),
                    0 if side == "3" else 1,
                ] = row.dG

    covered = ~np.isnan(watson_crick).any(axis=(2, 3))
    assert covered.any(), "the later fit contributed no phosphodiester cell"
    return watson_crick, mismatch


FINAL_WC, FINAL_MISMATCH = build_final_tables()

WC_ANCHOR = np.full((4, 4), -1, dtype=np.int8)
"""Row of the mismatch table for each Watson-Crick pair, -1 for a pair that is not one."""
for _pair, _row in zip(_WC_PAIRS, range(4)):
    WC_ANCHOR[DNA_BASES.index(_pair[0]), DNA_BASES.index(_pair[1])] = _row


def _mismatch_means(mismatch):
    """Mean mismatch cost per sugar pair, for a pair the fit has no specific cell for.

    A sugar pair the fit does not cover at all falls back to the mean over every cell.
    """
    overall = float(mismatch[np.isfinite(mismatch)].mean())
    means = np.full((3, 3), overall)
    for top in range(3):
        for bottom in range(3):
            cells = mismatch[top, bottom][np.isfinite(mismatch[top, bottom])]
            if cells.size:
                means[top, bottom] = cells.mean()
    return means


MISMATCH_MEAN = _mismatch_means(FINAL_MISMATCH)
"""Cost of a mismatch with no Watson-Crick neighbour to lean on, or one the fit misses."""


def build_step_energy():
    """Watson-Crick step energy for every sugar and base combination, kcal/mol.

    Indexed `[top1, top2, bottom1, bottom2, base, base]` over the four sugars the step spans.
    The later fit supplies the cells it has, which is the chemically uniform steps it was
    measured on, and the earlier one covers the rest through its class map, reverse
    complementing the dinucleotide where that fit pooled the two strand orientations.
    """
    energy = np.zeros((3, 3, 3, 3, 4, 4), dtype=np.float64)
    for top1 in range(3):
        for top2 in range(3):
            for bottom1 in range(3):
                for bottom2 in range(3):
                    cls, flip = CLASS_MAP[top1, top2, bottom1, bottom2]
                    for first in range(4):
                        for second in range(4):
                            value = np.nan
                            if top1 == top2 and bottom1 == bottom2:
                                value = FINAL_WC[top1, bottom1, first, second]
                            if np.isnan(value):
                                if flip:
                                    value = PARAMETERS[cls, COMPLEMENT[second], COMPLEMENT[first]]
                                else:
                                    value = PARAMETERS[cls, first, second]
                            energy[top1, top2, bottom1, bottom2, first, second] = value
    return energy


def build_mismatch_energy():
    """Cost of crossing one mismatched pair, kcal/mol.

    Indexed `[top sugar, bottom sugar, anchor base, anchor base, mismatch base, mismatch base]`.
    Only mismatches that have a Watson-Crick pair beside them are read from here; the fitted 3'
    and 5' entries are averaged when both exist, and a pair the fit has no cell for takes the
    mean cost for its sugar pair.
    """
    energy = np.zeros((3, 3, 4, 4, 4, 4), dtype=np.float64)
    for top in range(3):
        for bottom in range(3):
            fallback = MISMATCH_MEAN[top, bottom]
            for anchor_top in range(4):
                for anchor_bottom in range(4):
                    row = WC_ANCHOR[anchor_top, anchor_bottom]
                    for first in range(4):
                        for second in range(4):
                            value = fallback
                            if row >= 0:
                                three = FINAL_MISMATCH[top, bottom, row, first, second, 0]
                                five = FINAL_MISMATCH[top, bottom, row, first, second, 1]
                                if not np.isnan(three) and not np.isnan(five):
                                    value = 0.5 * (three + five)
                                elif not np.isnan(three):
                                    value = three
                                elif not np.isnan(five):
                                    value = five
                            energy[top, bottom, anchor_top, anchor_bottom, first, second] = value
    return energy


DINUCLEOTIDE_ENERGY = build_step_energy()
MISMATCH_ENERGY = build_mismatch_energy()


TANDEM_MISMATCH_ENERGY = MISMATCH_MEAN.copy()
"""Cost of a mismatch with no Watson-Crick pair beside it, kcal/mol, indexed `[top sugar, bottom sugar]`.

Two mismatches in a row are not two independent stacks, so neither the fitted anchored cells nor
their mean describes the pair. This currently holds the mean anchored cost for each sugar pair,
which is the finer-grained analogue of the single constant RIsearch1 charges every tandem
mismatch, and is the number to revisit. It takes no bases because nothing anchors it.
"""
