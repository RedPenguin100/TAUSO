"""Nearest-neighbour free-energy tables for DNA duplexes, loaded from CSV.

- PO/PO and PS/PO DNA duplexes: Wang SS, Xiong E, Bhadra S, Ellington AD,
  "Developing predictive hybridization models for phosphorothioate
  oligonucleotides using high-resolution melting" (PMC9116672). Reported as
  dG(50 C), dH, dS; the Gibbs free energy at 37 C is derived below.
- DNA/RNA hybrid: Sugimoto et al., Biochemistry 1995. dG37 (kcal/mol) plus the
  per-nearest-neighbour dH/dS behind it (used for melting-temperature calculations).

The tables live next to this module as CSV files (see weights/*.csv).
"""

from importlib.resources import files

import pandas as pd

from ...util import celsius_to_kelvin

_BODY_TEMPERATURE_K = celsius_to_kelvin(37.0)
_WEIGHTS = files("tauso.features.hybridization") / "weights"


def _load(filename):
    with (_WEIGHTS / filename).open() as handle:
        return pd.read_csv(handle, comment="#")


def _gibbs_37_weights(table):
    """Maps each dinucleotide to its dG37 (kcal/mol) from dH (kcal/mol) and dS (cal/mol/K)."""
    weights = {}
    for _, row in table.iterrows():
        weights[row["nucleotide"]] = row["H"] - _BODY_TEMPERATURE_K * (row["S"] / 1000.0)
    return weights


PO_PO_G37_WEIGHTS = _gibbs_37_weights(_load("po_po_dna_dna.csv"))  # Wang et al., phosphodiester DNA/DNA
PS_PO_G37_WEIGHTS = _gibbs_37_weights(_load("ps_po_dna_dna.csv"))  # Wang et al., phosphorothioate DNA/DNA

# Phosphorothioate backbone contribution at 37 C (kcal/mol), applied as a delta on top of the
# unmodified DNA/RNA hybrid baseline (DNA/DNA PO-vs-PS shift assumed transferable to the RNA-bound
# context). Positive = PS destabilises relative to the unmodified duplex.
PS_DELTA_DG37_WEIGHTS = {
    pair: PS_PO_G37_WEIGHTS[pair] - PO_PO_G37_WEIGHTS[pair]
    for pair in PO_PO_G37_WEIGHTS
    if pair in PS_PO_G37_WEIGHTS
}

# Sugimoto et al., DNA/RNA hybrid: dG37 (kcal/mol) and the dH/dS behind it.
_dna_rna = _load("dna_rna.csv")
DNA_RNA_DG37_WEIGHTS = dict(zip(_dna_rna["nucleotide"], _dna_rna["dG37"]))
DNA_RNA_DH_WEIGHTS = dict(zip(_dna_rna["nucleotide"], _dna_rna["dH"]))  # kcal/mol
DNA_RNA_DS_WEIGHTS = dict(zip(_dna_rna["nucleotide"], _dna_rna["dS"]))  # cal/(mol*K)
DNA_RNA_INIT = {"dH": 1.9, "dS": -3.9}  # Sugimoto helix-initiation factor (dH kcal/mol, dS cal/mol/K)
