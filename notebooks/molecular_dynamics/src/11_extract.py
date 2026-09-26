#!/usr/bin/env python
"""Step 11: measure the duplex geometry each trajectory holds.

Two modes, because the two halves cost very different things:

    11_extract.py --run <replicate dir>    analyse one trajectory -> features.csv beside it
    11_extract.py                          gather every features.csv -> out/features/perstep.csv

The per-run pass is a single cpptraj sweep and is what runs on the cluster, once per trajectory.
The gather is a concatenation and runs anywhere in seconds. Keeping them apart means a change to
how the measurements are pooled never costs 1,920 cpptraj passes again.

What is measured, and at which level it exists:

    pair    within one base pair: shear, stretch, stagger, buckle, propeller, opening, the
            hydrogen-bond count, the groove widths, and the C1'-C1' distance
    step    between consecutive pairs: shift, slide, rise, tilt, roll, twist, Zp, and the same
            geometry against a fitted local helical axis (hX, hY, hRise, hIncl, hTip, hTwist)
    res     per residue: sugar pucker and the backbone torsions, as circular statistics
    global  whole duplex: radius of gyration, and RMSD over all of it and over its core

The level is not a modelling choice -- a hydrogen bond belongs to a pair and Zp to the frame
between two pairs, and cpptraj writes them to separate files for that reason.

The first nanosecond of each production run is discarded. Block analysis over successive 1 ns
windows put it two to three times the between-window scatter away from the rest, in the direction
of the A-form starting geometry, with no residual drift after it.

nastruct needs to be told that the modXNA residue names are nucleotides; the map is generated from
the same tables step 3 built the libraries from, so it cannot drift from what was simulated.
"""
import argparse
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from consts import (BASES_BY_SUGAR, DUPLEX_LENGTH as L, FEATURES, MDIN, ROLES, SUGARS, SYSTEMS,
                    ensure_dirs)

PS_PER_FRAME = 2.0                  # how much time one saved frame stands for
EQUILIBRATION_FRAMES = 500          # the first ns, discarded; see the note above
MARGIN = 2                          # terminal pairs fray; the core masks below exclude them

PAIR_PARAMS = ("Shear", "Stretch", "Stagger", "Buckle", "Propeller", "Opening", "HB",
               "Major", "Minor")
STEP_PARAMS = ("Shift", "Slide", "Rise", "Tilt", "Roll", "Twist", "Zp")
HELIX_PARAMS = (("X-disp", "hX"), ("Y-disp", "hY"), ("Rise", "hRise"),
                ("Incl.", "hIncl"), ("Tip", "hTip"), ("Twist", "hTwist"))
COLUMNS = ["system", "replicate", "level", "index", "obs", "mean", "sd"]


def frame_interval_ps():
    """How much simulated time one saved frame covers, read from the production input."""
    text = (MDIN / "prod.in").read_text()
    step = float(re.search(r"\bdt\s*=\s*([0-9.]+)", text).group(1))
    every = int(re.search(r"\bntwx\s*=\s*([0-9]+)", text).group(1))
    return step * every


def check_frame_interval():
    """EQUILIBRATION_FRAMES counts frames, so it only means 1 ns while a frame is 2 ps."""
    interval = frame_interval_ps()
    if abs(interval - PS_PER_FRAME) > 1e-9:
        raise SystemExit(
            f"prod.in saves a frame every {interval:g} ps, not {PS_PER_FRAME:g}. "
            f"EQUILIBRATION_FRAMES = {EQUILIBRATION_FRAMES} would discard "
            f"{EQUILIBRATION_FRAMES * interval / 1000:g} ns rather than "
            f"{EQUILIBRATION_FRAMES * PS_PER_FRAME / 1000:g} ns.")


def resmap():
    """`resmap` arguments telling nastruct which base each modXNA residue carries."""
    return " ".join(
        f"resmap {prefix}{sugar}{letter}:{'C' if letter == '5' else letter}"
        for prefix, _, _ in ROLES.values()
        for sugar in SUGARS
        for letter in BASES_BY_SUGAR[sugar])


def cpptraj_input(topology, trajectory, out):
    """One sweep writing every measurement this step reads."""
    core = f":{MARGIN + 1}-{L - MARGIN},{L + MARGIN + 1}-{2 * L - MARGIN}"
    lines = [f"parm {topology}",
             f"trajin {trajectory} {EQUILIBRATION_FRAMES + 1} last",
             "autoimage",
             f"nastruct {resmap()} naout {out}/na.dat"]
    for i in range(1, 2 * L + 1):
        lines.append(f"pucker p{i} :{i}@C1' :{i}@C2' :{i}@C3' :{i}@C4' :{i}@O4' "
                     f"out {out}/pucker.dat range360")
    lines += [f"multidihedral MD alpha beta gamma delta epsilon zeta chin out {out}/tors.dat",
              f"radgyr rg :1-{2 * L} out {out}/rg.dat",
              f"rms rall first :1-{2 * L}&!@H= out {out}/rmsd.dat",
              f"rms rcore first {core}&!@H= out {out}/rmsd_core.dat"]
    for i in range(1, L + 1):
        lines.append(f"distance d{i} :{i}@C1' :{2 * L + 1 - i}@C1' out {out}/pairs.dat")
    return "\n".join(lines + ["go", "quit", ""])


def circular(degrees):
    """Mean angle and circular variance. A torsion near 0/360 has no meaningful linear mean."""
    radians = np.radians(degrees)
    cos, sin = np.cos(radians).mean(), np.sin(radians).mean()
    length = np.hypot(cos, sin)
    return np.degrees(np.arctan2(sin, cos)) % 360, 1 - length


def read(path):
    if not Path(path).exists():
        return None
    table = pd.read_csv(path, sep=r"\s+")
    table.columns = [c.lstrip("#") for c in table.columns]
    return table


def measure(work, out):
    """Every row this trajectory contributes, as (level, index, observable, mean, sd)."""
    rows = []

    table = read(out / "BP.na.dat")
    if table is not None:
        paired = table[table.Base1 + table.Base2 == 2 * L + 1]
        for position, group in paired.groupby("Base1"):
            for name in PAIR_PARAMS:
                if name in group:
                    rows.append(("pair", int(position), name,
                                 group[name].mean(), group[name].std()))

    for path, params in ((out / "BPstep.na.dat", [(p, p) for p in STEP_PARAMS]),
                         (out / "Helix.na.dat", list(HELIX_PARAMS))):
        table = read(path)
        if table is None:
            continue
        first = table.BP1.str.split("-").str[0].astype(int)
        second = table.BP1.str.split("-").str[1].astype(int)
        table = table.assign(i=first)[(first + second == 2 * L + 1).values]
        for position, group in table.groupby("i"):
            for column, name in params:
                if column in group:
                    rows.append(("step", int(position), name,
                                 group[column].mean(), group[column].std()))

    table = read(out / "pucker.dat")
    if table is not None:
        for column in table.columns:
            if column.startswith("p") and column[1:].isdigit():
                values = table[column].to_numpy(float)
                mean, variance = circular(values)
                rows.append(("res", int(column[1:]), "pucker", mean, variance))
                north = ((values < 36) | (values > 324)).mean()
                rows.append(("res", int(column[1:]), "north_frac", north, np.nan))

    table = read(out / "tors.dat")
    if table is not None:
        for column in table.columns:
            if ":" in column:
                name, index = column.split(":")
                mean, variance = circular(table[column].to_numpy(float))
                rows.append(("res", int(index), name, mean, variance))

    table = read(out / "pairs.dat")
    if table is not None:
        for column in table.columns:
            if column.startswith("d") and column[1:].isdigit():
                rows.append(("pair", int(column[1:]), "c1c1",
                             table[column].mean(), table[column].std()))

    for filename, name in (("rg.dat", "rg"), ("rmsd.dat", "rmsd"), ("rmsd_core.dat", "rmsd_core")):
        table = read(out / filename)
        if table is not None and table.shape[1] >= 2:
            values = table.iloc[:, 1].to_numpy(float)
            rows.append(("global", 0, name, values.mean(), values.std()))

    return rows


def analyse(work):
    """Measure one replicate directory. Returns the number of rows written."""
    system = work.parent
    topology, trajectory = system / "dry.prmtop", work / "md_nowat.nc"
    if not (topology.exists() and trajectory.exists()):
        raise SystemExit(f"missing {topology.name} or {trajectory.name} in {work}")

    scratch = work / ".extract"
    scratch.mkdir(exist_ok=True)
    result = subprocess.run(["cpptraj"], input=cpptraj_input(topology, trajectory, scratch),
                            text=True, capture_output=True, check=False)
    rows = measure(work, scratch)
    if not rows:
        raise SystemExit(f"cpptraj produced nothing for {work}\n{result.stdout[-2000:]}")

    frame = pd.DataFrame([(system.name, work.name) + r for r in rows], columns=COLUMNS)
    frame.to_csv(work / "features.csv", index=False)
    for leftover in scratch.iterdir():
        leftover.unlink()
    scratch.rmdir()
    return len(frame)


def gather():
    """Concatenate every per-run features.csv into the one table the tables step reads."""
    ensure_dirs()
    found = sorted(SYSTEMS.glob("*/*/features.csv"))
    if not found:
        raise SystemExit(f"no features.csv under {SYSTEMS} -- run with --run first")
    table = pd.concat([pd.read_csv(f) for f in found], ignore_index=True)
    out = FEATURES / "perstep.csv"
    table.to_csv(out, index=False)

    print(f"{len(table):,} measurements from {len(found):,} runs -> {out}")
    print(f"  duplexes   : {table.system.nunique()}")
    print(f"  replicates : {sorted(table.replicate.unique())}")
    for level, group in table.groupby("level"):
        print(f"  {level:7s} {group.obs.nunique():3d} observables x "
              f"{group['index'].nunique():3d} positions")
    return table


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run", type=Path, help="one replicate directory to analyse")
    args = parser.parse_args()

    if args.run:
        check_frame_interval()
        print(f"{args.run}: {analyse(args.run)} rows", flush=True)
    else:
        gather()


if __name__ == "__main__":
    sys.exit(main())
