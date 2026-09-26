"""Shared helpers for the pmemd stages (minimise, heat, equilibrate, production).

Kept in one place because every stage needs the same two awkward things: an environment with conda
stripped out, and a way to read Amber's output without picking up its trailing statistics. Both are
the kind of detail that is silently wrong rather than loudly broken, so there is one copy.
"""
import os
import re
import subprocess
import warnings

import parmed as pmd

from consts import AMBER_SH, MDIN, PMEMD, SOLVENT_RESIDUES

warnings.filterwarnings("ignore")

NAN_PATTERN = re.compile(r"=\s*NaN", re.I)
# Amber writes "<step> <energy> <rms> <gmax> ..." with no marker before the step number
ENERGY_LINE = re.compile(r"^\s*\d+\s+(-?\d+\.\d+E?[+-]?\d*)", re.M)
TEMP_LINE = re.compile(r"TEMP\(K\)\s*=\s*(-?\d+\.\d+)")
DENSITY_LINE = re.compile(r"Density\s*=\s*(-?\d+\.\d+)")

# Amber appends AVERAGES and RMS FLUCTUATIONS blocks after the last step, and both carry the same
# field names. Reading the last match in the file returns the RMS fluctuation, not the value.
AVERAGES_HEADER = "A V E R A G E S"


def pmemd_environment():
    """os.environ with conda stripped. pmemd inherits conda's libstdc++ otherwise and dies
    silently right after printing its CUDA banner -- no error, no NaN, just a missing restart."""
    env = os.environ.copy()
    for key in ("LD_LIBRARY_PATH", "PATH"):
        env[key] = os.pathsep.join(
            p for p in env.get(key, "").split(os.pathsep) if p and "conda" not in p.lower())
    for key in ("CONDA_PREFIX", "CONDA_DEFAULT_ENV", "PYTHONHOME", "PYTHONPATH"):
        env.pop(key, None)
    if AMBER_SH.exists():
        env["AMBERHOME"] = str(AMBER_SH.parent)
        env["LD_LIBRARY_PATH"] = os.pathsep.join(
            filter(None, [str(AMBER_SH.parent / "lib"), env.get("LD_LIBRARY_PATH")]))
    return env


def count_solute_residues(prmtop):
    """Leading non-solvent residues; restraint masks are stamped per system, never hardcoded."""
    topology = pmd.load_file(str(prmtop))
    count = 0
    for residue in topology.residues:
        if residue.name in SOLVENT_RESIDUES:
            break
        count += 1
    return count


def last_value(out_file, pattern):
    """The last value at a real step, ignoring Amber's trailing statistics blocks."""
    steps_only = out_file.read_text().split(AVERAGES_HEADER)[0]
    values = pattern.findall(steps_only)
    return float(values[-1]) if values else None


def write_input(stage, work, n_solute, seed=None):
    """Stamp a control file for one system: its own restraint mask, and its replicate's seed.

    The mask has to name the residues this system actually has, and minimisation has no seed to
    stamp, so both are substituted here rather than written into config/md.
    """
    text = (MDIN / f"{stage}.in").read_text().replace("SOLUTE", str(n_solute))
    if seed is not None:
        text = text.replace("SEED", str(seed))
    (work / f"{stage}.in").write_text(text)


def run_stage(work, stage, topology, coords, env, reference=None, trajectory=None):
    """One pmemd call inside a replicate directory. Returns (ok, detail)."""
    restart = work / f"{stage}.rst7"
    if restart.exists():
        return True, "already done"

    command = [str(PMEMD), "-O", "-i", f"{stage}.in", "-p", topology, "-c", coords,
               "-r", restart.name, "-o", f"{stage}.out"]
    if reference:
        command += ["-ref", reference]
    if trajectory:
        command += ["-x", trajectory]

    subprocess.run(command, cwd=work, env=env, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL, check=False)

    out = work / f"{stage}.out"
    if out.exists() and NAN_PATTERN.search(out.read_text()):
        return False, f"NaN during {stage} -- bad ion draw, re-solvate this system"
    if not restart.exists():
        return False, f"{stage} wrote no restart -- check {stage}.out"
    return True, "ok"
