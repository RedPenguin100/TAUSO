#!/usr/bin/env python
"""Step 10: the 4 ns production run -- the only trajectory that is analysed.

Everything before this was preparation. Here the duplex runs free at 310 K and 1 bar with no
restraints, and a frame is written every 2 ps.

    ntr=0          no restraints at all. Any bias here would show up directly in the observable,
                   which is the fluctuation of a base pair.
    ntwx=1000      a frame every 2 ps -> 2000 frames from 4 ns. That is the sampling the C1'-C1'
                   standard deviation is computed over.
    ioutfm=1       binary NetCDF rather than ASCII
    iwrap=1        keep molecules inside the periodic box, so nothing appears to fly apart

Run length is not settled by a convergence test. 4 ns with five replicates is chosen to sample
more than one trajectory can, not because a pilot showed 4 ns is where a base-pair fluctuation
converges -- that pilot, a handful of duplexes to 50 ns with the weights recomputed from
2/5/10/20/50 ns windows, has not been run. Worth knowing before these are presented as final.

Storage. A solvated frame is every water in the box, so one run's prod.nc is ~900 MB and the whole
campaign is ~1.7 TB -- more than the scratch holds. Stripping the solvent leaves ~26 MB, which is
all the analysis reads. With KEEP_SOLVATED False the solvated trajectory is deleted once the
stripped one is confirmed to hold every frame it should; the deletion is irreversible, and what is
lost is the ability to re-examine solvent structure later.
"""
import re
import subprocess
import sys

from consts import MDIN, PMEMD, REPLICATE_SEEDS, SYSTEMS
from mdrun import (TEMP_LINE, count_solute_residues, last_value, pmemd_environment, run_stage,
                   write_input)

KEEP_SOLVATED = False            # keep prod.nc as well as md_nowat.nc; see the note above
SOLVENT_MASK = ":WAT,Na+,Cl-"
FRAMES_READ = re.compile(r"Read (\d+) frames")


def expected_frames():
    """How many frames prod.in writes: one every ntwx steps over nstlim."""
    text = (MDIN / "prod.in").read_text()
    nstlim = int(re.search(r"nstlim\s*=\s*(\d+)", text).group(1))
    ntwx = int(re.search(r"ntwx\s*=\s*(\d+)", text).group(1))
    return nstlim // ntwx


def strip_water(system, work):
    """Drop solvent from the trajectory. Returns the number of frames the stripped copy holds.

    cpptraj needs the solvated topology to read prod.nc and `strip` applied after `trajin`; the
    dry topology is written separately because the analysis reads the two together and a stripped
    trajectory against a solvated topology silently reads garbage.
    """
    dry_top = system / "dry.prmtop"
    if not dry_top.exists():
        subprocess.run(["cpptraj"], input=(
            f"parm {system / 'sys.prmtop'}\n"
            f"parmstrip {SOLVENT_MASK}\n"
            f"parmwrite out {dry_top}\nrun\nquit\n"),
            text=True, capture_output=True, check=False)

    out = work / "md_nowat.nc"
    if not out.exists():
        result = subprocess.run(["cpptraj"], input=(
            f"parm {system / 'sys.prmtop'}\n"
            f"trajin {work / 'prod.nc'}\nautoimage\nstrip {SOLVENT_MASK}\n"
            f"trajout {out}\nrun\nquit\n"),
            text=True, capture_output=True, check=False)
        if not out.exists():
            return 0
        found = FRAMES_READ.search(result.stdout)
        return int(found.group(1)) if found else 0

    # already stripped on an earlier pass; count what is there rather than trusting the file
    result = subprocess.run(["cpptraj"], input=(
        f"parm {dry_top}\ntrajin {out}\nrun\nquit\n"), text=True, capture_output=True, check=False)
    found = FRAMES_READ.search(result.stdout)
    return int(found.group(1)) if found else 0


def produce_one(system, replicate, seed, env, wanted):
    label = f"{system.name}/{replicate}"
    work = system / replicate
    if not (work / "equil.rst7").exists():
        return label, False, "not equilibrated -- run 09_equilibrate.py first"

    if not (work / "prod.rst7").exists():
        write_input("prod", work, count_solute_residues(system / "sys.prmtop"), seed)
        ok, detail = run_stage(work, "prod", "../sys.prmtop", "equil.rst7", env,
                               trajectory="prod.nc")
        if not ok:
            return label, False, detail

    frames = strip_water(system, work)
    if frames != wanted:
        return label, False, f"stripped trajectory has {frames} frames, expected {wanted}"

    solvated = work / "prod.nc"
    if not KEEP_SOLVATED and solvated.exists():
        solvated.unlink()

    return label, True, last_value(work / "prod.out", TEMP_LINE)


def main(limit=None):
    if not PMEMD.exists():
        raise SystemExit(f"pmemd not found at {PMEMD}")
    if not (MDIN / "prod.in").exists():
        raise SystemExit(f"no prod.in in {MDIN}")

    systems = sorted(p for p in SYSTEMS.iterdir() if p.is_dir())[:limit or None]
    if not systems:
        raise SystemExit(f"no systems in {SYSTEMS} -- run 06_solvate.py first")

    env = pmemd_environment()
    wanted = expected_frames()
    total = len(systems) * len(REPLICATE_SEEDS)
    print(f"production: {len(systems)} systems x {len(REPLICATE_SEEDS)} replicates = {total} runs")
    print(f"  {wanted} frames each, ~20 min per run on one GPU -- this is the expensive step")
    print(f"  solvated trajectory is {'kept' if KEEP_SOLVATED else 'deleted once stripped'}")

    results = []
    for system in systems:
        for replicate, seed in REPLICATE_SEEDS.items():
            results.append(produce_one(system, replicate, seed, env, wanted))
            label, ok, detail = results[-1]
            print(f"    {label}: {'ok' if ok else detail}", flush=True)

    failures = [(label, detail) for label, ok, detail in results if not ok]
    temps = [d for _, ok, d in results if ok and d is not None]

    print(f"\n  produced : {len(results) - len(failures)}/{total}")
    if temps:
        print(f"  temperature : {min(temps):.1f} to {max(temps):.1f} K")
    if failures:
        nans = [f for f in failures if "NaN" in str(f[1])]
        print(f"\n  FAILURES ({len(failures)}, of which {len(nans)} are NaN):")
        for label, reason in failures[:10]:
            print(f"    {label}: {reason}")
        raise SystemExit(1)
    print("\n  all trajectories complete")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else None)
