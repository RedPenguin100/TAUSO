"""Paths and shared constants for the nearest-neighbour pipeline.

Everything resolves from this file's own location, so the project can be moved or cloned anywhere
without editing a path. Import it rather than hardcoding directories:

    from consts import PDBS, RESIDUES, SYSTEMS

Layout:
    src/     this code
    config/  what to run: the design (which sequences, which chemistries) and the MD control files
    lib/     external dependencies (the pinned modXNA clone)
    out/     everything generated: residues, structures, systems, trajectories
"""
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

SRC = ROOT / "src"
CONFIG = ROOT / "config"
LIB = ROOT / "lib"
OUT = ROOT / "out"

# --- what to run ---------------------------------------------------------------------------------
MDIN = CONFIG / "md"                      # min1.in, min2.in, heat.in, ...
DESIGN = CONFIG / "design" / "design.json"   # written by 01 and 02, read by everything after
MODXNA = LIB / "modxna"                   # pinned clone; commit is what every result was built with
MODXNA_SH = MODXNA / "modxna.sh"
FRCMOD = MODXNA / "dat" / "frcmod.modxna"

# --- generated ----------------------------------------------------------------------------------
RESIDUES = OUT / "residues"               # 03: the 48 .lib templates
PDBS = OUT / "pdbs"                       # 04: starting coordinates
VALIDATE = OUT / "validate"               # 05: dry tleap builds and their checks
SYSTEMS = OUT / "systems"                 # 06+: solvated systems, then all MD output
FEATURES = OUT / "features"               # 11: the per-step fluctuations the weights are fitted to

# --- engines ------------------------------------------------------------------------------------
# pmemd lives outside conda on purpose; see pmemd_environment in mdrun.py about libstdc++
PMEMD = Path("/home/michael/molecular_dynamics_new/pmemd24/bin/pmemd.cuda")
AMBER_SH = PMEMD.parent.parent / "amber.sh"

# --- chemistry ----------------------------------------------------------------------------------
DUPLEX_LENGTH = 16                        # base pairs
COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}

SUGARS = {                                # letter -> (modXNA fragment, label)
    "D": ("DC2", "2'-deoxy"),
    "M": ("MOE", "2'-MOE"),
    "E": ("CET", "cEt"),
    "R": ("RC3", "ribose"),
}
RIBOSE_LIKE = "R"                        # sugars whose backbone and bases are the RNA ones
BASES = {"A": "DAA", "G": "DGG", "T": "DTT", "C": "DCC", "5": "M5C"}
BASES_RNA = {"A": "RAA", "G": "RGG", "U": "RUU", "C": "RCC"}
# Cytosine is 5-methyl on every synthetic strand, which is what a synthesised gapmer carries.
# A ribose strand stands for the RNA target and keeps plain cytosine, since a transcript is not
# methylated at C5. The rule follows the sugar at each position, not the strand as a whole.
BASES_BY_SUGAR = {"D": "AGT5", "M": "AGT5", "E": "AGT5", "R": "AGUC"}
ROLES = {                                 # role -> (residue-name prefix, backbone fragment, flag)
    "internal": ("N", "DPO", ""),
    "cap5":     ("Q", "5PO", "--5cap"),
    "cap3":     ("Z", "DPO", "--3cap"),
}


def backbone_fragment(role, sugar):
    """Ribose carries its own phosphate fragment; the 5' cap has no phosphate whatever the sugar."""
    return "RPO" if sugar in RIBOSE_LIKE and role != "cap5" else ROLES[role][1]


def base_fragment(sugar, letter):
    """RNA spells thymine as uracil, and its bases are separate modXNA fragments."""
    return (BASES_RNA if sugar in RIBOSE_LIKE else BASES)[letter]


def base_letter(sugar, letter):
    """Spell one base in the alphabet its own sugar uses.

    Sequences are written in the plain DNA alphabet throughout, so this is where a position picks
    up what its chemistry actually carries: a ribose position reads T as U, and a synthetic
    position reads C as 5, the 5-methyl cytosine a synthesised strand is made with.
    """
    if sugar in RIBOSE_LIKE:
        return "U" if letter == "T" else letter
    return "5" if letter == "C" else letter


TARGET_CHARGE = {"N": -1.000, "Q": -0.320, "Z": -0.680}

SOLVENT_RESIDUES = {"WAT", "HOH", "Na+", "Cl-", "K+"}

# every tleap build sources exactly these, dry or solvated: the structure step 5 checks has to be
# the structure step 10 integrates, and a force field that differed between them would not show up
# in any check either one runs
FORCE_FIELD = ("leaprc.DNA.OL21", "leaprc.RNA.OL3", "leaprc.water.opc")

# Replicates differ in one thing only: the random seed used to assign starting velocities at
# heating. Solvation and minimisation are shared, so they live above the replicate directories.
REPLICATE_SEEDS = {"rep1": 11111, "rep2": 22222, "rep3": 33333,
                   "rep4": 44444, "rep5": 55555}
MINUTES_PER_RUN = 20             # 4 ns on the 5090; used only for the estimate below


def describe_design(rows):
    """Print what a design contains: duplexes per cell, and what that costs to run."""
    print(f"{len(rows)} duplexes -> {DESIGN}")
    replicates = len(REPLICATE_SEEDS)
    for (family, state), n in sorted(Counter((r["family"], r["state"]) for r in rows).items()):
        print(f"  {family:13s} {state:16s} {n:4d} duplexes -> {n * replicates:5d} runs")
    runs = len(rows) * replicates
    print(f"\n  {runs} runs at {replicates} replicates, "
          f"{runs * MINUTES_PER_RUN / 60:.0f} GPU-hours at ~{MINUTES_PER_RUN} min each")


def reverse_complement(sequence):
    return "".join(COMPLEMENT[b] for b in reversed(sequence))


def ensure_dirs():
    for path in (CONFIG, MDIN, DESIGN.parent, LIB, OUT, RESIDUES, PDBS, VALIDATE,
                 SYSTEMS, FEATURES):
        path.mkdir(parents=True, exist_ok=True)
