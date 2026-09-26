#!/usr/bin/env python
"""Step 6: put each duplex in water and salt.

A duplex in vacuum is unphysical -- a -30 charge with nothing screening it. Four tleap commands fix
that, and the order matters:

    addions m Na+ 0            neutralise. "0" means "as many as needed", not zero. Run BEFORE
                               solvation so counterions land on the phosphates by electrostatics
                               rather than being dropped in bulk water. PME assumes a neutral cell.
    solvateoct m OPCBOX 10.0   truncated octahedron, ~30% fewer waters than a cube for the same
                               minimum image distance. The buffer must exceed the nonbonded cutoff,
                               or the duplex interacts with its own periodic image.
    addionsrand m Na+ <n>      salt, to 150 mM. `addionsrand` takes a COUNT, not a concentration,
    addionsrand m Cl- <n>      and the count that gives 150 mM depends on how much water the box
                               holds -- which is not known until the box exists. So the system is
                               solvated twice: once to learn its water count, then again with the
                               count that puts it at 150 mM.

addionsrand places ions at RANDOM water sites, so every solvation is a different draw and this step
is not reproducible from the script alone -- sys.rst7 is the only record of what was built. Roughly
3% of draws put an ion somewhere dynamically unstable: fine at the minimum, then the system explodes
on heating. That is handled by re-solvating later, not here.
"""
import json
import re
import shutil
import subprocess
import warnings
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

from consts import DESIGN, PDBS, RESIDUES, SYSTEMS, ensure_dirs
from residue_library import leap_preamble

warnings.filterwarnings("ignore")

BUFFER_ANGSTROM = 10.0
SALT_MOLAR = 0.15                # physiological ionic strength
WATER_MOLAR = 55.5               # molarity of pure water, to convert a count into a concentration
WORKERS = 8
PROGRESS_EVERY = 32
ADDED_RESIDUES = re.compile(r"Added (\d+) residues")

# first pass: solvate only, to find out how much water this particular box holds
COUNT_INPUT = """\
{preamble}
m = loadpdb {pdb}
addions m Na+ 0
solvateoct m OPCBOX {buffer}
quit
"""

LEAP_INPUT = """\
{preamble}
m = loadpdb {pdb}
addions m Na+ 0
solvateoct m OPCBOX {buffer}
addionsrand m Na+ {salt}
addionsrand m Cl- {salt}
saveamberparm m sys.prmtop sys.rst7
quit
"""


def ion_pairs(water):
    """How many NaCl pairs put a box of `water` waters at SALT_MOLAR.

    `addionsrand` replaces water with ions rather than adding to it, so each pair removes two
    waters from the very total the concentration is measured against. Solving for that leaves
    n = C*W / (55.5 + 2C) rather than the naive C*W/55.5.
    """
    return round(SALT_MOLAR * water / (WATER_MOLAR + 2 * SALT_MOLAR))


def count_water(out, pdb):
    """Solvate once without salt and read the water count off tleap's own report."""
    (out / "count.in").write_text(COUNT_INPUT.format(
        preamble=leap_preamble(RESIDUES), pdb=pdb, buffer=BUFFER_ANGSTROM))
    with open(out / "count.log", "w") as log:
        subprocess.run(["tleap", "-f", "count.in"], cwd=out,
                       stdout=log, stderr=subprocess.STDOUT, check=False)
    found = ADDED_RESIDUES.findall((out / "count.log").read_text())
    return int(found[-1]) if found else 0


def solvate_one(pdb):
    out = SYSTEMS / pdb.stem
    out.mkdir(parents=True, exist_ok=True)
    if (out / "sys.prmtop").exists():
        return pdb.stem, out, True                 # already done; the script is resumable

    water = count_water(out, pdb)
    if not water:
        return pdb.stem, out, False
    (out / "solvate.in").write_text(LEAP_INPUT.format(
        preamble=leap_preamble(RESIDUES), pdb=pdb,
        buffer=BUFFER_ANGSTROM, salt=ion_pairs(water)))

    with open(out / "solvate.log", "w") as log:
        subprocess.run(["tleap", "-f", "solvate.in"], cwd=out,
                       stdout=log, stderr=subprocess.STDOUT, check=False)
    return pdb.stem, out, (out / "sys.prmtop").exists()


def inspect(out):
    """What tleap actually built, read straight off the prmtop's flag sections.

    Only five numbers are wanted, so the topology is never assembled into objects: building 37,000
    atoms and 9,000 residues per system costs a couple of seconds each, and there are 384 of them.
    Charges are stored pre-multiplied by 18.2223, so summing first and converting once is both
    faster and closer than converting each of them.
    """
    import parmed as pmd

    data = pmd.amber.AmberFormat(str(out / "sys.prmtop")).parm_data
    counts = Counter(data["RESIDUE_LABEL"])
    water = counts["WAT"]
    return dict(atoms=data["POINTERS"][0], charge=sum(data["CHARGE"]) / 18.2223,
                water=water, na=counts["Na+"], cl=counts["Cl-"],
                box=round(data["BOX_DIMENSIONS"][1], 1),
                salt=round(counts["Cl-"] / (water / 55.5), 2) if water else 0.0)


def main():
    if shutil.which("tleap") is None:
        raise SystemExit("tleap not on PATH -- conda activate aso_tools")

    known = {r["name"] for r in json.loads(DESIGN.read_text())}
    pdbs = [p for p in sorted(PDBS.glob("*.pdb")) if p.stem in known]
    if not pdbs:
        raise SystemExit(f"no PDBs in {PDBS} -- run 04_build_pdbs.py first")
    ensure_dirs()

    print(f"solvating {len(pdbs)} structures "
          f"({BUFFER_ANGSTROM} A buffer, NaCl to {SALT_MOLAR * 1000:.0f} mM)")
    results = []
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        for done, result in enumerate(pool.map(solvate_one, pdbs), start=1):
            results.append(result)
            if done % PROGRESS_EVERY == 0 or done == len(pdbs):
                print(f"  {done}/{len(pdbs)}", flush=True)

    failures, summary = [], []
    for name, out, ok in results:
        if not ok:
            failures.append((name, "tleap produced no topology"))
            continue
        info = inspect(out)
        summary.append(info)
        if abs(info["charge"]) > 1e-3:              # PME needs a neutral cell
            failures.append((name, f"system not neutral: {info['charge']:+.3f}"))

    print(f"\n  solvated : {len(summary)}/{len(pdbs)}")
    if summary:
        def span(key):
            return min(x[key] for x in summary), max(x[key] for x in summary)

        print(f"  atoms     : {span('atoms')[0]} to {span('atoms')[1]}")
        print(f"  waters    : {span('water')[0]} to {span('water')[1]}")
        print(f"  Na+ / Cl- : {span('na')[0]}-{span('na')[1]} / {span('cl')[0]}-{span('cl')[1]}")
        print(f"  box edge  : {span('box')[0]} to {span('box')[1]} A")
        print(f"  salt      : ~{span('salt')[0]} to {span('salt')[1]} M")
        print(f"  charge    : {max(abs(x['charge']) for x in summary):.1e} (all neutral)")

    if failures:
        print(f"\n  FAILURES ({len(failures)}):")
        for name, why in failures[:10]:
            print(f"    {name}: {why}")
        raise SystemExit(1)
    print("\n  all systems neutral and solvated")


if __name__ == "__main__":
    main()
