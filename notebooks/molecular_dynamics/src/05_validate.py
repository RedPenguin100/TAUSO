#!/usr/bin/env python
"""Step 5: build each structure with tleap and prove it is correct, before spending any compute.

tleap does the real assembly -- it resolves the residue names written in step 4 against the
libraries from step 3, grows the atoms the PDB never contained (the 2'-MOE arm, the cEt bridge,
the 5-methyl on every cytosine, every hydrogen) and applies the force field. The checks then read
properties off the object MD will actually integrate.

    tleap errors      an unknown residue means the PDB and the libraries disagree on a name;
                      a missing parameter means the force field cannot describe something
    solute molecules  must be exactly 2. A fragmented backbone -- strands broken into loose
                      residues -- runs to completion and NEVER NaNs, yielding a clean-looking
                      trajectory of free-floating nucleotides. This count is the only place it
                      gets caught, and it silently invalidated an earlier campaign.
    net charge        must be an integer, or the electrostatics are quietly wrong in every frame
    delta torsion     ~83 deg is C3'-endo (A-form), ~145 is C2'-endo; confirms fd_helix geometry
                      actually landed rather than wc_helix
    atom counts       every residue must have exactly the atoms its library template defines.
                      Each strand comes from its own fibre model, and a residue that arrives short
                      of a base is rebuilt by tleap on top of its neighbour -- which passes every
                      other check here while quietly destroying the pairing.
"""
import json
import shutil
import subprocess
import warnings
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from consts import DESIGN, PDBS, RESIDUES, VALIDATE, ensure_dirs
from residue_library import leap_preamble, read_lib

warnings.filterwarnings("ignore")

LEAP_INPUT = """\
{preamble}
m = loadpdb {pdb}
saveamberparm m vac.prmtop vac.rst7
quit
"""

TLEAP_ERRORS = ("fatal", "unknown residue", "could not find bond",
                "could not find angle", "could not find dihedral")
SUGAR_TORSION = ("C5'", "C4'", "C3'", "O3'")
WORKERS = 8


def build_one(pdb):
    """Run tleap on one structure. Returns (name, workdir, errors found in the log)."""
    work = VALIDATE / pdb.stem
    work.mkdir(parents=True, exist_ok=True)

    (work / "tleap.in").write_text(
        LEAP_INPUT.format(preamble=leap_preamble(RESIDUES), pdb=pdb))

    with open(work / "tleap.log", "w") as log:
        subprocess.run(["tleap", "-f", "tleap.in"], cwd=work,
                       stdout=log, stderr=subprocess.STDOUT, check=False)

    log_text = (work / "tleap.log").read_text().lower()
    return pdb.stem, work, [e for e in TLEAP_ERRORS if e in log_text]


def count_molecules(top):
    """Connected components of the bond graph. Two strands should be two molecules."""
    graph = [[] for _ in top.atoms]
    for bond in top.bonds:
        graph[bond.atom1.idx].append(bond.atom2.idx)
        graph[bond.atom2.idx].append(bond.atom1.idx)

    seen, molecules = set(), 0
    for start in range(len(top.atoms)):
        if start in seen:
            continue
        molecules += 1
        seen.add(start)
        stack = [start]
        while stack:
            for neighbour in graph[stack.pop()]:
                if neighbour not in seen:
                    seen.add(neighbour)
                    stack.append(neighbour)
    return molecules


def sugar_torsions(top):
    coordinates = top.coordinates        # a property that rebuilds an array on every access
    angles = []
    for residue in top.residues:
        atoms = {a.name: a.idx for a in residue.atoms}
        if not all(name in atoms for name in SUGAR_TORSION):
            continue
        p = [coordinates[atoms[name]] for name in SUGAR_TORSION]
        b0, b1, b2 = p[1] - p[0], p[2] - p[1], p[3] - p[2]
        n1, n2 = np.cross(b0, b1), np.cross(b1, b2)
        angles.append(np.degrees(np.arctan2(
            np.dot(np.cross(n1, n2), b1 / np.linalg.norm(b1)), np.dot(n1, n2))) % 360)
    return angles


def template_atom_counts():
    """How many atoms each residue name should have, read off the libraries built in step 3."""
    return {path.stem: len(read_lib(path)[0]) for path in RESIDUES.glob("*.lib")}


def inspect(work, expected):
    import parmed as pmd

    top = pmd.load_file(str(work / "vac.prmtop"), str(work / "vac.rst7"))
    angles = sugar_torsions(top)
    wrong = [f"{r.name}{r.idx + 1} has {len(r.atoms)}, expects {expected[r.name]}"
             for r in top.residues if r.name in expected and len(r.atoms) != expected[r.name]]
    return dict(molecules=count_molecules(top), residues=len(top.residues),
                charge=sum(a.charge for a in top.atoms), wrong_atoms=wrong,
                torsion_min=min(angles), torsion_max=max(angles))


def main():
    if shutil.which("tleap") is None:
        raise SystemExit("tleap not on PATH -- conda activate aso_tools")

    known = {r["name"] for r in json.loads(DESIGN.read_text())}
    built = sorted(PDBS.glob("*.pdb"))
    pdbs = [p for p in built if p.stem in known]
    stale = len(built) - len(pdbs)
    if stale:
        print(f"ignoring {stale} PDBs that are not in the current design")
    if not pdbs:
        raise SystemExit(f"no PDBs in {PDBS} -- run 04_build_pdbs.py first")
    ensure_dirs()

    expected = template_atom_counts()
    print(f"validating {len(pdbs)} structures against {len(expected)} residue templates")
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        results = list(pool.map(build_one, pdbs))

    failures, summary = [], []
    for name, work, errors in results:
        if errors:
            failures.append((name, ", ".join(errors)))
            continue
        if not (work / "vac.prmtop").exists():
            failures.append((name, "tleap produced no topology"))
            continue

        info = inspect(work, expected)
        summary.append(info)
        if info["molecules"] != 2:
            failures.append((name, f"{info['molecules']} molecules, expected 2 -- FRAGMENTED"))
        elif abs(info["charge"] - round(info["charge"])) > 1e-3:
            failures.append((name, f"non-integer charge {info['charge']:.3f}"))
        elif info["wrong_atoms"]:
            failures.append((name, "; ".join(info["wrong_atoms"][:2])))

    print(f"\n  built and parsed : {len(summary)}/{len(pdbs)}")
    if summary:
        def tally(values):
            return dict(Counter(values))

        print(f"  solute molecules : {tally(x['molecules'] for x in summary)}"
              f"   (2 = both strands intact)")
        print(f"  residues         : {tally(x['residues'] for x in summary)}")
        print(f"  net charge       : {tally(round(x['charge']) for x in summary)}")
        print(f"  delta torsion    : {min(x['torsion_min'] for x in summary):.1f} to "
              f"{max(x['torsion_max'] for x in summary):.1f} deg   (~83 = C3'-endo, A-form)")
        print(f"  atom counts      : "
              f"{sum(1 for x in summary if not x['wrong_atoms'])}/{len(summary)} match templates")

    if failures:
        print(f"\n  FAILURES ({len(failures)}):")
        for name, reason in failures[:10]:
            print(f"    {name}: {reason}")
        raise SystemExit(1)
    print("\n  all checks passed")


if __name__ == "__main__":
    main()
