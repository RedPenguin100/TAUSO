#!/usr/bin/env python
"""Step 13: decide which observables become features.

An observable earns a place only if a sixteen-cell dinucleotide table predicts it, that prediction
beats what base composition alone gives, and it is not a restatement of something already kept.
This script computes all three and writes the table the paper reports.

    predictability   variance of the observable explained by its dinucleotide cell means,
                     1 - within-cell variance / total variance
    composition      the same quantity against the G+C count of the dinucleotide, three levels.
                     The control that separates a structural feature from a GC proxy
    redundancy       correlation between two observables across the sixteen cell values. A pair
                     counts as redundant only if it exceeds REDUNDANT in EVERY chemistry: two
                     descriptors can coincide in one state and diverge in another, and only the
                     former is a genuine duplicate

Predictability is fitted to the mean over an observable's five replicates rather than to single
runs, so that measurement noise is not charged against it; the reproducibility of those averages
is reported alongside.

Every figure is computed per chemistry state and reported at its worst, since a feature is applied
to all of them. Junction states are screened separately from the uniform arms: a junction step is
the boundary between two chemistries and belongs to neither.
"""
import json
import sys
from itertools import combinations

import numpy as np
import pandas as pd

from consts import DESIGN, DUPLEX_LENGTH as L, FEATURES, ensure_dirs

MARGIN = 2                  # terminal steps fray and are excluded, as everywhere else
MIN_PREDICTABILITY = 0.40   # in every chemistry
REDUNDANT = 0.90            # |r| across the sixteen cells
UNIFORM = {"D:R": "deoxy gap", "M:R": "2'-MOE wing", "E:R": "cEt wing"}


def explained(values, key):
    """Fraction of variance the cell means account for. Zero when the key says nothing."""
    key, values = np.asarray(key), np.asarray(values, float)
    if len(values) < 50 or values.var() == 0:
        return np.nan
    means = np.array([values[key == k].mean() for k in key])
    return 1 - ((values - means) ** 2).mean() / values.var()


def annotate(design):
    """Per-step rows with the dinucleotide, its G+C count, and which chemistry the step sits in."""
    table = pd.read_csv(FEATURES / "perstep.csv")
    table = table[(table.level == "step") &
                  (table["index"] > MARGIN) & (table["index"] < L - MARGIN)].copy()

    def sugars(name):
        spec = design[name]["sugarA"]
        return spec * L if len(spec) == 1 else spec

    table["dinucleotide"] = [design[s]["seq"][i - 1:i + 1]
                             for s, i in zip(table.system, table["index"])]
    table["gc"] = table.dinucleotide.map(lambda d: sum(b in "GC" for b in d))
    # a step whose two positions carry different sugars is a junction, and gets its own state
    table["state"] = [
        f"{sugars(s)[i-1]}{sugars(s)[i]}|RR" if sugars(s)[i - 1] != sugars(s)[i]
        else design[s]["state"]
        for s, i in zip(table.system, table["index"])]
    return table


def screen(table, states):
    """Predictability, composition and reproducibility per observable, at its worst state."""
    rows = []
    for observable, group in table.groupby("obs"):
        record = {"observable": observable}
        worst = {}
        for state in states:
            here = group[group.state == state]
            if here.empty:
                continue
            noise = here.groupby(["system", "index"])["mean"].var().mean()
            cell = here.groupby(["system", "index"]).agg(
                value=("mean", "mean"), dinucleotide=("dinucleotide", "first"), gc=("gc", "first"))
            nn = explained(cell.value, cell.dinucleotide)
            gc = explained(cell.value, cell.gc)
            gc = 0.0 if np.isnan(gc) else gc
            record[f"{state}_nn"] = round(nn, 2)
            record[f"{state}_gc"] = round(gc, 2)
            worst[state] = (nn, gc, 1 - (noise / 5) / cell.value.var())
        if not worst:
            continue
        record["nn_worst"] = round(min(v[0] for v in worst.values()), 2)
        record["gain_worst"] = round(min(v[0] - v[1] for v in worst.values()), 2)
        record["reproducibility"] = round(min(v[2] for v in worst.values()), 2)
        rows.append(record)
    return pd.DataFrame(rows)


def correlations(table, states):
    """Every observable pair's correlation in every chemistry, as a tidy table.

    Kept because the redundancy decision turns on the weakest of these, not the strongest, and a
    reader cannot check that from the verdict alone.
    """
    rows = []
    for state in states:
        wide = table[table.state == state].groupby(["obs", "dinucleotide"])["mean"].mean().unstack(0)
        if wide.shape[0] < 16:
            continue
        for a, b in combinations(wide.columns, 2):
            rows.append(dict(state=state, a=a, b=b, r=round(abs(wide[a].corr(wide[b])), 3)))
    return pd.DataFrame(rows)


def redundancy(table, states):
    """For each observable, its strongest partner among those it duplicates in every chemistry."""
    seen = {}
    for state in states:
        wide = table[table.state == state].groupby(["obs", "dinucleotide"])["mean"].mean().unstack(0)
        if wide.shape[0] < 16:
            continue
        for a, b in combinations(wide.columns, 2):
            seen.setdefault((a, b), []).append(abs(wide[a].corr(wide[b])))

    worst = {}
    for (a, b), correlations in seen.items():
        if len(correlations) < len(states) or min(correlations) <= REDUNDANT:
            continue                       # coincides in some states but not all: not a duplicate
        weakest = min(correlations)
        for x, y in ((a, b), (b, a)):
            if weakest > worst.get(x, (0, None))[0]:
                worst[x] = (weakest, y)
    return worst


def verdict(row, partner):
    """Why an observable is kept or dropped, in the order the reasons are applied."""
    r, other = partner.get(row["observable"], (0, None))
    if other:
        return f"redundant with {other} in every chemistry (|r| >= {r:.2f})"
    if row["nn_worst"] < MIN_PREDICTABILITY:
        return f"predictability {row['nn_worst']:.2f} below {MIN_PREDICTABILITY}"
    return "retained"


def report(name, table, states):
    result = screen(table, states)
    partner = redundancy(table, states)
    # a redundant pair keeps one member; the helical-frame descriptor is the one retained
    keep = {a for a in partner if a.startswith("h")}
    result["verdict"] = [
        "retained" if row["observable"] in keep else verdict(row, partner)
        for row in result.to_dict("records")]
    result["closest"] = [partner.get(o, (None, None))[1] for o in result.observable]
    result["closest_r_min"] = [partner.get(o, (None, None))[0] for o in result.observable]
    result = result.sort_values("gain_worst", ascending=False)

    print(f"\n{name}\n")
    head = "".join(f"{s:>14s}" for s in states)
    print(f"{'observable':11s}{head}{'worst':>8s}{'gain':>7s}{'repro':>7s}   verdict")
    print(f"{'':11s}" + "".join(f"{'NN':>7s}{'GC':>7s}" for _ in states))
    for row in result.to_dict("records"):
        cells = "".join(f"{row[f'{s}_nn']:7.2f}{row[f'{s}_gc']:7.2f}" for s in states)
        print(f"{row['observable']:11s}{cells}{row['nn_worst']:8.2f}{row['gain_worst']:7.2f}"
              f"{row['reproducibility']:7.2f}   {row['verdict']}")
    return result


def main():
    ensure_dirs()
    design = {r["name"]: r for r in json.loads(DESIGN.read_text())}
    table = annotate(design)

    uniform = report("Uniform arms", table, list(UNIFORM))
    junctions = sorted({s for s in table.state.unique() if "|" in s})
    junction = report("Junction states", table, junctions)

    out = FEATURES / "screen.csv"
    pd.concat([uniform.assign(arm="uniform"), junction.assign(arm="junction")]).to_csv(
        out, index=False)
    pairs = FEATURES / "correlations.csv"
    pd.concat([correlations(table, list(UNIFORM)),
               correlations(table, junctions)]).to_csv(pairs, index=False)
    retained = uniform[uniform.verdict == "retained"].observable.tolist()
    print(f"\n{len(retained)} retained: {', '.join(sorted(retained))}")
    print(f"-> {out}")
    print(f"-> {pairs}")

    wide = pd.read_csv(FEATURES / "correlations.csv")
    wide = wide[wide.state.isin(UNIFORM)]
    span = wide.groupby(["a", "b"]).r.agg(["min", "max"])
    split = span[(span["max"] > REDUNDANT) & (span["min"] <= REDUNDANT)]
    if len(split):
        print(f"\npairs redundant in some chemistries but not all -- both kept:")
        for (a, b), row in split.sort_values("max", ascending=False).head(6).iterrows():
            print(f"  {a:7s} ~ {b:7s}  |r| {row['min']:.2f} to {row['max']:.2f}")


if __name__ == "__main__":
    sys.exit(main())
