#!/usr/bin/env python
"""Step 4: build starting coordinates for every duplex.

Geometry comes from nab's fd_helix canonical A-form fibre model. This choice is load-bearing: the
obvious alternative, wc_helix, lays down a rigid B-form duplex regardless of what it is asked for --
including for RNA, which is never B-form. A run of this length cannot cross the sugar-pucker
barrier, so the trajectory reports back whatever pucker it was built with. Months of runs once
concluded 2'-MOE sits C2'-endo; that was the builder, not the chemistry. Do not swap this back.

fd_helix returns a duplex whose second strand is the exact complement of the first. Every duplex
in this campaign is Watson-Crick, so that is what is wanted -- but the assembly does not assume
it: strand A comes from the model of `seq`, strand B from the model of `partner`, superposed onto
the place the first left for it. Whatever `partner` says is what gets built, and no base is ever
rebuilt from scratch. The backbone of an A-form fibre strand does not depend on its sequence, so
the superposition is exact rather than approximate.

Every sugar in the fibre model carries a 2'-OH, and that oxygen is the handle for every chemistry:

    2'-deoxy   delete it
    2'-MOE     keep it -- the methoxyethyl arm attaches there
    cEt        rename it O6' -- in cEt that same oxygen bridges C2' to C4'
    ribose     keep it -- the fibre model is already RNA, so nothing is touched at all

Only the anchor is placed; tleap grows the rest of the arm or bridge from the library templates.
Hydrogens are dropped for the same reason.
"""
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from consts import DESIGN, DUPLEX_LENGTH as L, PDBS, ROLES, base_letter, ensure_dirs

NAB_SOURCE = r"""
// read "outfile sequence" pairs on stdin, write one canonical A-form duplex per line
string name, seq;
molecule m;

while ( scanf("%s %s", name, seq) == 2 ) {
    m = fd_helix("arna", seq, "rna");
    putpdb(name, m, "-wwpdb");
}
"""
# atoms shared by every residue whatever its base, so the frame two strands are matched on
BACKBONE = ("P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'")


def build_nab(workdir):
    """nab is a compiled language: build the binary once, then stream every duplex through it."""
    source, binary = workdir / "fd.nab", workdir / "fd"
    source.write_text(NAB_SOURCE)
    subprocess.run(["nab",
                    "--compiler", os.environ.get("NAB_CC", "gcc"),
                    "--linker", os.environ.get("NAB_FC", "gfortran"),
                    str(source), "-o", str(binary)],
                   cwd=workdir, check=True, capture_output=True, text=True)
    return binary


def read_atoms(path):
    """Heavy atoms of a PDB as (residue index, atom name, xyz, template line)."""
    out = []
    with open(path) as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            name = line[12:16].strip()
            if name.startswith("H"):
                continue                                  # tleap rebuilds hydrogens
            out.append((int(line[22:26]), name,
                        np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                        line))
    return out


def superpose(mobile, target):
    """Rotation and translation putting the mobile points on the target ones."""
    mc, tc = mobile.mean(axis=0), target.mean(axis=0)
    u, _, vt = np.linalg.svd((mobile - mc).T @ (target - tc))
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rotation = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return rotation, tc - rotation @ mc


def strand_frame(atoms, first, last):
    """Backbone coordinates of residues first..last, keyed so two strands can be matched up."""
    return {(res - first, name): xyz for res, name, xyz, _ in atoms
            if first <= res <= last and name in BACKBONE}


def residue_name(position, sugar, base):
    role = "cap5" if position == 1 else "cap3" if position == L else "internal"
    return f"{ROLES[role][0]}{sugar}{base_letter(sugar, base)}"


def per_position(spec):
    """A sugar given once applies to every position; given per position it is used as written."""
    return spec * L if len(spec) == 1 else spec


def write_line(line, serial, resname, resseq, xyz):
    return (f"ATOM  {serial:>5d} {line[12:16]}{line[16]}{resname:>3}{line[20:22]}{resseq:>4d}"
            f"{line[26]}   {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{line[54:]}")


def build_duplex(row, geometry):
    """Assemble one duplex: strand A as built, strand B taken from its own model and moved here."""
    a_atoms = read_atoms(geometry[row["seq"]])
    b_atoms = read_atoms(geometry[row["partner"]])

    # strand B belongs where strand 2 of A's own model sits; both are A-form strands of the same
    # length, so matching them on the backbone is exact.
    here = strand_frame(a_atoms, L + 1, 2 * L)
    there = strand_frame(b_atoms, 1, L)
    shared = sorted(set(here) & set(there))
    mobile = np.array([there[k] for k in shared])
    target = np.array([here[k] for k in shared])
    rotation, shift = superpose(mobile, target)
    fit = np.sqrt(np.mean(np.sum((mobile @ rotation.T + shift - target) ** 2, axis=1)))

    sugars = (per_position(row["sugarA"]), per_position(row["sugarB"]))
    bases = (row["seq"], row["partner"])
    lines, serial = [], 0
    for strand, atoms in enumerate((a_atoms, b_atoms)):
        for res, name, xyz, line in atoms:
            if not 1 <= res <= L:
                continue                                  # each model contributes one strand
            sugar = sugars[strand][res - 1]
            if name == "O2'":
                if sugar == "D":
                    continue                              # 2'-deoxy has no 2' oxygen
                if sugar == "E":
                    name, line = "O6'", line[:12] + " O6'" + line[16:]  # in cEt it bridges to C4'
            if strand:
                xyz = rotation @ xyz + shift
            serial += 1
            lines.append(write_line(
                line, serial,
                residue_name(res, sugar, bases[strand][res - 1]),
                res + strand * L, xyz))
    return "".join(lines) + "TER\nEND\n", fit


def main():
    if shutil.which("nab") is None:
        raise SystemExit("nab not on PATH -- activate the environment that has it:\n"
                         "    conda activate aso_build")
    if not DESIGN.exists():
        raise SystemExit(f"no design at {DESIGN} -- run 01_design_wc.py first")

    ensure_dirs()
    # duplexes already simulated keep their structures wherever they were run; only what is
    # still outstanding needs building here.
    rows = [r for r in json.loads(DESIGN.read_text()) if r.get("status", "todo") == "todo"]
    strands = sorted({s for row in rows for s in (row["seq"], row["partner"])})

    failed, fits = [], []
    with tempfile.TemporaryDirectory(prefix="fdhelix_") as tmp:
        workdir = Path(tmp)
        binary = build_nab(workdir)
        # fd_helix builds an A-form fibre over ACGU, so T is fed as U. Thymine's methyl and
        # cytosine's 5-methyl are both atoms tleap grows from the residue library, so neither
        # has to reach nab; only the backbone and the ring positions come from the fibre.
        paths = {s: workdir / f"geom{i:04d}.pdb" for i, s in enumerate(strands)}
        subprocess.run([str(binary)], cwd=workdir, text=True, capture_output=True, check=True,
                       input="".join(f"{p} {s.replace('T', 'U').lower()}\n"
                                     for s, p in paths.items()))

        for row in rows:
            if not all(paths[s].exists() for s in (row["seq"], row["partner"])):
                failed.append((row["name"], "fd_helix produced nothing"))
                continue
            text, fit = build_duplex(row, paths)
            fits.append(fit)
            (PDBS / f"{row['name']}.pdb").write_text(text)

    print(f"built {len(rows) - len(failed)}/{len(rows)} outstanding PDBs "
          f"from {len(strands)} unique strands -> {PDBS}")
    if fits:
        print(f"  strand B superposition: {min(fits):.4f} to {max(fits):.4f} A "
              f"(exact, both strands are the same fibre model)")
    for name, reason in failed[:8]:
        print(f"  FAILED {name}: {reason}")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
