#!/usr/bin/env python
"""Step 9: let the box reach its real density at 310 K and 1 bar.

Heating ran at constant volume, so the box is still whatever size tleap made it -- the system is at
temperature but not at the right density. Equilibration switches on the barostat and lets the volume
settle, which is the point of the stage.

    ntb=2, ntp=1   constant pressure; the box is now allowed to change size
    ntx=5, irest=1 continue from heating WITH its velocities, rather than reassigning them
    ntt=3          Langevin thermostat, the same one production uses -- switching thermostat
                   between equilibration and production would mean equilibrating the wrong dynamics
    restraint 0.5  down from 25 during heating: the duplex is nudged, not held, so it can relax
                   into the new box while still being discouraged from drifting or fraying
    ntwx=0         no trajectory. Nothing here is analysed; this is 500 ps of settling.
"""
import sys

from consts import MDIN, PMEMD, REPLICATE_SEEDS, SYSTEMS
from mdrun import (DENSITY_LINE, TEMP_LINE, count_solute_residues, last_value,
                   pmemd_environment, run_stage, write_input)

PROGRESS_EVERY = 100


def equilibrate_one(system, replicate, seed, env):
    label = f"{system.name}/{replicate}"
    work = system / replicate
    if not (work / "heat.rst7").exists():
        return label, False, "not heated -- run 08_heat.py first"
    if (work / "equil.rst7").exists():
        return label, True, (last_value(work / "equil.out", TEMP_LINE),
                             last_value(work / "equil.out", DENSITY_LINE))

    write_input("equil", work, count_solute_residues(system / "sys.prmtop"), seed)
    ok, detail = run_stage(work, "equil", "../sys.prmtop", "heat.rst7", env,
                           reference="heat.rst7")
    if not ok:
        return label, False, detail
    return label, True, (last_value(work / "equil.out", TEMP_LINE),
                         last_value(work / "equil.out", DENSITY_LINE))


def main(limit=None):
    if not PMEMD.exists():
        raise SystemExit(f"pmemd not found at {PMEMD}")
    if not (MDIN / "equil.in").exists():
        raise SystemExit(f"no equil.in in {MDIN}")

    systems = sorted(p for p in SYSTEMS.iterdir() if p.is_dir())[:limit or None]
    if not systems:
        raise SystemExit(f"no systems in {SYSTEMS} -- run 06_solvate.py first")

    env = pmemd_environment()
    total = len(systems) * len(REPLICATE_SEEDS)
    print(f"equilibrating {len(systems)} systems x {len(REPLICATE_SEEDS)} replicates "
          f"= {total} runs")

    results = []
    for system in systems:
        for replicate, seed in REPLICATE_SEEDS.items():
            results.append(equilibrate_one(system, replicate, seed, env))
            if len(results) % PROGRESS_EVERY == 0 or len(results) == total:
                print(f"  {len(results)}/{total}", flush=True)

    failures = [(label, detail) for label, ok, detail in results if not ok]
    good = [detail for _, ok, detail in results if ok and detail[0] is not None]

    print(f"\n  equilibrated : {len(results) - len(failures)}/{total}")
    if good:
        temps = [t for t, _ in good]
        densities = [d for _, d in good if d is not None]
        print(f"  temperature : {min(temps):.1f} to {max(temps):.1f} K   (target 310)")
        if densities:
            print(f"  density     : {min(densities):.3f} to {max(densities):.3f} g/mL   "
                  f"(water is ~1.0)")
    if failures:
        print(f"\n  FAILURES ({len(failures)}):")
        for label, reason in failures[:10]:
            print(f"    {label}: {reason}")
        raise SystemExit(1)
    print("\n  all replicates equilibrated")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else None)
