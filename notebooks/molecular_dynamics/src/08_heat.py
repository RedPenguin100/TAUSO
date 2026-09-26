#!/usr/bin/env python
"""Step 8: bring each minimised system from rest up to 310 K.

After minimisation the system sits at 0 K with no velocities at all. Heating assigns them from a
Maxwell distribution and ramps the thermostat 100 -> 310 K over the first 10 ps, then holds for 90.

    ntb=1, ntp=0   constant volume. Density is not settled yet; that is equilibration's job, and
                   letting the box breathe while the system is still cold gives a bad density.
    ntc=2, ntf=2   SHAKE on bonds to hydrogen, which is what makes a 2 fs timestep stable.
    ntr=1, wt=25   the solute is held firmly so water and ions equilibrate around a fixed duplex
                   before it is allowed to move.

This is the stage that fails. Roughly 3% of runs NaN here, and the cause is upstream: addionsrand
placed an ion somewhere dynamically unstable -- fine at the minimum, then the electrostatics
collapse on heating. It is not fixable by precision, thermostat or ramp rate; the fix is a different
random draw, meaning re-solvate that system and redo it. This script says so rather than letting a
broken system flow into production.

Replicates begin here. They differ in the velocity seed and nothing else, so solvation and
minimisation are shared and only the dynamics live in per-replicate directories.
"""
import sys

from consts import MDIN, PMEMD, REPLICATE_SEEDS, SYSTEMS
from mdrun import (TEMP_LINE, count_solute_residues, last_value, pmemd_environment, run_stage,
                   write_input)

PROGRESS_EVERY = 100


def heat_one(system, replicate, seed, env):
    """Heat one system for one replicate. Returns (label, ok, detail)."""
    label = f"{system.name}/{replicate}"
    if not (system / "min2.rst7").exists():
        return label, False, "not minimised -- run 07_minimize.py first"

    work = system / replicate
    work.mkdir(parents=True, exist_ok=True)
    if (work / "heat.rst7").exists():
        return label, True, last_value(work / "heat.out", TEMP_LINE)

    write_input("heat", work, count_solute_residues(system / "sys.prmtop"), seed)
    ok, detail = run_stage(work, "heat", "../sys.prmtop", "../min2.rst7", env,
                           reference="../min2.rst7", trajectory="heat.nc")
    if not ok:
        return label, False, detail
    return label, True, last_value(work / "heat.out", TEMP_LINE)


def main(limit=None):
    if not PMEMD.exists():
        raise SystemExit(f"pmemd not found at {PMEMD}")
    if not (MDIN / "heat.in").exists():
        raise SystemExit(f"no heat.in in {MDIN}")

    systems = sorted(p for p in SYSTEMS.iterdir() if p.is_dir())[:limit or None]
    if not systems:
        raise SystemExit(f"no systems in {SYSTEMS} -- run 06_solvate.py first")

    env = pmemd_environment()
    total = len(systems) * len(REPLICATE_SEEDS)
    print(f"heating {len(systems)} systems x {len(REPLICATE_SEEDS)} replicates = {total} runs")

    results = []
    for system in systems:
        for replicate, seed in REPLICATE_SEEDS.items():
            results.append(heat_one(system, replicate, seed, env))
            if len(results) % PROGRESS_EVERY == 0 or len(results) == total:
                print(f"  {len(results)}/{total}", flush=True)

    failures = [(label, detail) for label, ok, detail in results if not ok]
    temperatures = [detail for _, ok, detail in results if ok and detail is not None]

    print(f"\n  heated : {len(results) - len(failures)}/{total}")
    if temperatures:
        print(f"  final temperature : {min(temperatures):.1f} to {max(temperatures):.1f} K "
              f"(target 310)")
    if failures:
        nans = [f for f in failures if "NaN" in str(f[1])]
        print(f"\n  FAILURES ({len(failures)}, of which {len(nans)} are NaN):")
        for label, reason in failures[:10]:
            print(f"    {label}: {reason}")
        raise SystemExit(1)
    print("\n  all replicates at temperature")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else None)
