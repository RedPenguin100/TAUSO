#!/usr/bin/env python
"""Step 3: build the residue definitions the force field does not ship.

Amber knows deoxyribose and ribose. It has never heard of a 2'-methoxyethyl arm or a 2',4'-bridged
sugar, so each modified residue has to be defined: atom names and types, partial charges,
connectivity, internal geometry. modXNA assembles one from three published fragments --

    backbone  +  sugar  +  base   ->  residue
    DPO          MOE       DAA        NMA, an internal 2'-MOE adenosine

modXNA ships only the fragments; the residues are generated here, so this script plus the pinned
clone reproduces them exactly.

Residues are role-specific because the ends of a strand differ chemically: a 5' terminus has no
phosphate, a 3' terminus carries one plus a free hydroxyl. Different atom count, different charge.

Phosphodiester rather than phosphorothioate. Cytosine is 5-methyl on every synthetic strand and
plain on ribose, which is what a synthesised gapmer bound to a transcript carries -- so the
methylated residues `ND5`, `NM5`, `NE5` and their caps are built here alongside the rest.
"""
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor

from consts import (BASES_BY_SUGAR, MODXNA_SH, RESIDUES, ROLES, SUGARS, TARGET_CHARGE,
                    backbone_fragment, base_fragment, ensure_dirs)
from residue_library import read_lib

REQUIRED_TOOLS = ("cpptraj", "tleap", "sander")     # modXNA drives all three
JOBS_PER_RUN = 8


def make_jobs():
    """Every residue to build: (name, fragment spec, cap flag)."""
    return [(f"{prefix}{sugar}{base}",
             f"{backbone_fragment(role, sugar)} {fragment} {base_fragment(sugar, base)}", cap)
            for role, (prefix, _, cap) in ROLES.items()
            for sugar, (fragment, _) in SUGARS.items()
            for base in BASES_BY_SUGAR[sugar]]


def build_one(job):
    name, spec, cap = job
    workdir = RESIDUES / f".work_{name}"
    shutil.rmtree(workdir, ignore_errors=True)
    workdir.mkdir(parents=True)

    (workdir / f"in_{name}.txt").write_text(spec + "\n")
    command = ["bash", str(MODXNA_SH), "-i", f"in_{name}.txt", "-m", name] + ([cap] if cap else [])

    with open(RESIDUES / f"log_{name}.txt", "w") as log:
        subprocess.run(command, cwd=workdir, stdout=log, stderr=subprocess.STDOUT,
                       timeout=1800, check=False)

    built = workdir / f"{name}.lib"
    ok = built.exists()
    if ok:
        shutil.move(str(built), RESIDUES / f"{name}.lib")
    shutil.rmtree(workdir, ignore_errors=True)
    return name, ok


def check_charges(names):
    """A residue off by even 0.03 gives the duplex a non-integer net charge, and the electrostatics
    are then quietly wrong in every frame of every run. Cheapest place to catch it."""
    charges = {n: read_lib(RESIDUES / f"{n}.lib")[1] for n in names}
    errors = [(n, round(c, 4)) for n, c in charges.items()
              if abs(c - TARGET_CHARGE[n[0]]) > 1e-3]
    if errors:
        print(f"charge check FAILED: {errors}")
        return False
    print("charge check: all exact")
    return True


def print_atom_counts():
    labels = [label for _, label in SUGARS.values()]
    print("\n" + " " * 10 + "".join(f"{label:>10s}" for label in labels) + "   atoms (adenine)")
    for role, (prefix, _, _) in ROLES.items():
        counts = [len(read_lib(RESIDUES / f"{prefix}{sugar}A.lib")[0]) for sugar in SUGARS]
        print(f"{role:10s}" + "".join(f"{n:>10d}" for n in counts))


def main():
    if not MODXNA_SH.exists():
        raise SystemExit(f"modXNA not found at {MODXNA_SH} -- clone it into lib/ first")
    missing = [t for t in REQUIRED_TOOLS if shutil.which(t) is None]
    if missing:
        raise SystemExit(f"missing on PATH: {', '.join(missing)}\n"
                         f"modXNA needs the Amber tools:\n    conda activate aso_tools")

    ensure_dirs()
    jobs = make_jobs()
    per_sugar = " + ".join(f"{s}:{len(b)}" for s, b in BASES_BY_SUGAR.items())
    print(f"{len(jobs)} residues to have ({per_sugar} bases, x {len(ROLES)} roles)")

    todo = [j for j in jobs if not (RESIDUES / f"{j[0]}.lib").exists()]
    print(f"  {len(jobs) - len(todo)} already built, building {len(todo)}")
    with ThreadPoolExecutor(max_workers=JOBS_PER_RUN) as pool:
        results = list(pool.map(build_one, todo))

    failed = [name for name, ok in results if not ok]
    if failed:
        raise SystemExit(f"FAILED to build: {failed}")

    names = [j[0] for j in jobs]
    print(f"built {len(results)}, have {len(names)}")
    if not check_charges(names):
        raise SystemExit(1)
    print_atom_counts()
    print(f"\nwrote {len(names)} libraries to {RESIDUES}")


if __name__ == "__main__":
    main()
