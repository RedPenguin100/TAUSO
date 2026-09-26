#!/usr/bin/env python
"""Step 7: minimise each solvated system in two stages.

solvateoct drops water onto a grid, so the built system has bad contacts everywhere -- waters
overlapping the solute, ions too close to phosphates -- plus strain inside the duplex from atoms
tleap grew off templates. Heating that directly turns a 0.4 A overlap into kinetic energy and the
run explodes.

    min1   solute restrained, solvent free   water and ions find their places first, so the solute
                                             cannot deform around a badly placed water and keep it
    min2   everything free                   the solute relieves its own strain

min2 is load-bearing, not polish: every 2'-MOE and cEt structure is built with its 5'-cap arm
clashing, and min2 is what swings the arm clear. Shortening it reintroduces a heat-time explosion
that looks like a thermostat problem and is not one.

Minimisation is not dynamics -- no temperature, no time, just downhill until forces are small. It
is also shared by every replicate, since replicates differ only in the velocities assigned at
heating, so it is done once per system and lives above the replicate directories.
"""
import sys

from consts import MDIN, PMEMD, SYSTEMS
from mdrun import ENERGY_LINE, count_solute_residues, pmemd_environment, run_stage, write_input

PROGRESS_EVERY = 25


def minimise(system, env):
    """Both minimisation stages for one system. Returns (name, ok, detail)."""
    name = system.name
    topology = system / "sys.prmtop"
    if not topology.exists():
        return name, False, "no solvated topology -- run 06_solvate.py first"

    # the restraint mask has to match this system, so it is stamped in rather than hardcoded
    n_solute = count_solute_residues(topology)
    for stage in ("min1", "min2"):
        write_input(stage, system, n_solute, seed=None)

    # min1: solute restrained, so water and ions relax around a duplex held in place
    ok, detail = run_stage(system, "min1", "sys.prmtop", "sys.rst7", env, reference="sys.rst7")
    if not ok:
        return name, False, detail

    # min2: everything free, so the solute relieves its own strain
    ok, detail = run_stage(system, "min2", "sys.prmtop", "min1.rst7", env)
    if not ok:
        return name, False, detail

    energies = ENERGY_LINE.findall((system / "min2.out").read_text())
    if not energies:
        return name, False, "could not read a final energy from min2.out"
    return name, True, (n_solute, float(energies[-1]))


def main(limit=None):
    if not PMEMD.exists():
        raise SystemExit(f"pmemd not found at {PMEMD}")
    if not (MDIN / "min1.in").exists():
        raise SystemExit(f"no min1.in in {MDIN}")

    systems = sorted(p for p in SYSTEMS.iterdir() if p.is_dir())[:limit or None]
    if not systems:
        raise SystemExit(f"no systems in {SYSTEMS} -- run 06_solvate.py first")

    env = pmemd_environment()
    print(f"minimising {len(systems)} systems on {PMEMD.name}")

    results = []
    for done, system in enumerate(systems, start=1):
        results.append(minimise(system, env))
        if done % PROGRESS_EVERY == 0 or done == len(systems):
            print(f"  {done}/{len(systems)}", flush=True)

    failures = [(name, detail) for name, ok, detail in results if not ok]
    successes = [detail for _, ok, detail in results if ok]

    print(f"\n  minimised : {len(successes)}/{len(systems)}")
    if successes:
        print(f"  solute residues restrained in min1 : {sorted({s[0] for s in successes})}")
        energies = [s[1] for s in successes]
        print(f"  final energy : {min(energies):.0f} to {max(energies):.0f} kcal/mol")

    if failures:
        print(f"\n  FAILURES ({len(failures)}):")
        for name, reason in failures[:10]:
            print(f"    {name}: {reason}")
        raise SystemExit(1)
    print("\n  all systems minimised")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else None)
