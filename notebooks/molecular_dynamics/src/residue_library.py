"""The residue libraries step 3 builds: reading one, and loading them all into tleap.

The libraries step 3 builds are the reference for what a residue is: which atoms it has and what
they are charged. Step 3 checks the charges as it builds, step 5 checks that tleap grew exactly
the atoms the template defines, and both read the format here.

A .lib is a flat text file of `!entry` sections. The one that matters is `.unit.atoms`, whose
lines run

    "name" "type" typex resx flags seq elmnt charge

with the charge LAST -- the field before it is the atomic number, which is easy to grab by mistake.
"""
import re
from functools import cache

from consts import FORCE_FIELD, FRCMOD

ATOM_LINE = re.compile(r'\s*"(\S+)"\s+"(\S+)"\s+\d+\s+\d+\s+\d+\s+\d+\s+(-?\d+)\s+(-?[\d.]+)')


def read_lib(path):
    """Atom names and total charge of a .lib unit."""
    names, charge, reading = [], 0.0, False
    with open(path) as handle:
        for line in handle:
            if line.startswith("!entry") and ".unit.atoms " in line:
                reading = True
                continue
            if not reading:
                continue
            if line.startswith("!"):
                break
            match = ATOM_LINE.match(line)
            if match:
                names.append(match.group(1))
                charge += float(match.group(4))
    return names, charge


@cache
def leap_preamble(residues):
    """The tleap header every build shares: force field, modXNA parameters, residue libraries."""
    return "\n".join([*(f"source {name}" for name in FORCE_FIELD),
                       f"loadamberparams {FRCMOD}",
                       *(f"loadoff {path}" for path in sorted(residues.glob("*.lib")))])
