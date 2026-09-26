#!/usr/bin/env python
"""Step 12: pool the per-position measurements into nearest-neighbour tables.

One table per structural level. A cell is a (chemistry state, sequence key) pair, and its value is
the mean over every duplex, replicate and position that falls in it.

    step   keyed by dinucleotide          a base-pair step is a property of the two pairs it
                                          spans, and the data agrees: the dinucleotide explains
                                          0.89 of Roll in the gap against 0.14 for a single base
    pair   keyed by base, and by triplet   the pair alone is NOT enough -- propeller reaches 0.45
                                          by base and 0.78 by triplet -- so both keys are written,
                                          and the triplet is the one to prefer where it has support
    res    keyed by base                   pucker and the backbone torsions belong to one residue

Junction steps get their own cells. A step whose two positions carry different sugars is not
described by either uniform state, and a gapmer contains two of them.

Terminal positions are excluded: the two pairs at each end fray, and report the fraying rather
than the stacking.

Error bars resample DUPLEXES, not replicates. Five replicates of one duplex differ only in their
starting velocities, so their spread measures sampling noise; what limits a cell is how few
distinct sequence contexts it rests on.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from consts import DESIGN, FEATURES, ensure_dirs

PERSTEP = FEATURES / "perstep.csv"
OUT = FEATURES / "tables"

L = 16
MARGIN = 2
BOOTSTRAP = 400
MIN_DUPLEXES = 3
rng = np.random.default_rng(7)


def sugars_at(row):
    """Per-position sugar for strand A; strand B is ribose throughout."""
    spec = row["sugarA"]
    return spec * L if len(spec) == 1 else spec


def cells_for(row):
    """(level, index) -> (group, sub) for every position this duplex contributes."""
    seq, partner, sa = row["seq"], row["partner"], sugars_at(row)
    out = {}
    for i in range(MARGIN + 1, L - MARGIN + 1):            # pair positions
        out[("pair", i)] = (f"{sa[i-1]}:R", seq[i-1])
        if 1 < i < L:
            out[("pair3", i)] = (f"{sa[i-1]}:R", seq[i-2:i+1])
    for i in range(MARGIN + 1, L - MARGIN):                # step i spans pairs i, i+1
        a, b = sa[i-1], sa[i]
        group = f"{a}:R" if a == b else f"{a}{b}|RR"
        out[("step", i)] = (group, seq[i-1:i+1])
    for i in range(MARGIN + 1, L - MARGIN + 1):            # residues, both strands
        out[("res", i)] = (f"{sa[i-1]}:R", seq[i-1])
        out[("res", L + i)] = ("R:R", partner[i-1])
    return out


def bootstrap(frame):
    """Spread of the cell mean under resampling of its duplexes."""
    duplexes = frame.system.unique()
    if len(duplexes) < MIN_DUPLEXES:
        return np.nan
    by_duplex = {d: g["mean"].to_numpy() for d, g in frame.groupby("system")}
    draws = [np.concatenate([by_duplex[d] for d in rng.choice(duplexes, len(duplexes))]).mean()
             for _ in range(BOOTSTRAP)]
    return float(np.std(draws))


def main():
    ensure_dirs()
    OUT.mkdir(parents=True, exist_ok=True)
    design = {r["name"]: r for r in json.loads(DESIGN.read_text())}
    d = pd.read_csv(PERSTEP)
    d = d[d.level != "global"]

    keys = {name: cells_for(row) for name, row in design.items()}
    assign = [keys.get(s, {}).get((lv, ix)) for s, lv, ix in
              zip(d.system, d.level, d["index"])]
    d["group"] = [a[0] if a else None for a in assign]
    d["sub"] = [a[1] if a else None for a in assign]

    # the triplet key reads the same pair rows under a wider context
    trip = d[d.level == "pair"].copy()
    tassign = [keys.get(s, {}).get(("pair3", ix)) for s, ix in zip(trip.system, trip["index"])]
    trip["group"] = [a[0] if a else None for a in tassign]
    trip["sub"] = [a[1] if a else None for a in tassign]
    trip["level"] = "pair3"
    d = pd.concat([d, trip], ignore_index=True)
    d = d[d.group.notna() & d["sub"].notna() & d["mean"].notna()]

    rows = []
    for (level, obs, group, sub), g in d.groupby(["level", "obs", "group", "sub"]):
        for stat, column in (("mean", "mean"), ("sd", "sd")):
            values = g[column].dropna()
            if values.empty:
                continue
            rows.append(dict(level=level, obs=obs, stat=stat, cell=f"{group}@{sub}",
                             group=group, sub=sub, value=values.mean(),
                             err=bootstrap(g) if stat == "mean" else np.nan,
                             n_positions=len(g), n_duplexes=g.system.nunique(),
                             n_replicates=g.replicate.nunique()))
    table = pd.DataFrame(rows)

    for level, name in (("step", "step.csv"), ("pair", "pair.csv"),
                        ("pair3", "pair_triplet.csv"), ("res", "res.csv")):
        part = table[table.level == level].drop(columns="level")
        part.to_csv(OUT / name, index=False)
        cells = part[part.stat == "mean"]
        print(f"{name:26s} {len(part):6,} rows  {cells.cell.nunique():4d} cells  "
              f"{cells.obs.nunique():3d} observables  "
              f"{cells.n_duplexes.min():3d}-{cells.n_duplexes.max():3d} duplexes/cell")
    print(f"\nstates: {sorted(table.group.unique())}")


if __name__ == "__main__":
    main()
