"""Nearest-neighbour parameters for phosphodiester (PO) and phosphorothioate (PS) DNA duplexes.

Source: Hassan et al., Nucleic Acids Research 2022, PMC9116672.
https://pmc.ncbi.nlm.nih.gov/articles/PMC9116672/

Table 1 (PO/PO) and Table 3 (PS/PO) are transcribed verbatim below as
``nucleotide, dG(50C), dH, dS`` with dG/dH in kcal/mol and dS in cal/(mol*K),
adjusted to 1 M NaCl and 50 C. ``E`` rows are the paper's terminal-end
parameters; they are retained for completeness and are simply never matched when
summing over the A/C/G/U dinucleotides of a sequence.

The phosphorothioate backbone effect used as a feature is the per-dinucleotide
difference ``ps_po - po_po`` evaluated at 37 C (positive = PS destabilises the
duplex relative to the unmodified phosphodiester).
"""

from io import StringIO

import pandas as pd

from ...util import celsius_to_kelvin

_BODY_TEMPERATURE_K = celsius_to_kelvin(37.0)

# Table 1 — phosphodiester / phosphodiester (PO/PO)
po_po_table = """
AA,-0.83,-8.10,-22.5
UU,-0.83,-8.10,-22.5
AU,-0.56,-5.53,-15.4
UA,-0.58,-6.40,-18.0
CA,-0.95,-6.89,-18.4
UG,-0.95,-6.89,-18.4
GU,-0.94,-7.12,-19.1
AC,-0.94,-7.12,-19.1
CU,-0.94,-7.51,-20.3
AG,-0.94,-7.51,-20.3
GA,-0.88,-6.51,-17.4
UC,-0.88,-6.51,-17.4
CG,-1.62,-10.81,-28.5
GC,-1.76,-12.68,-33.8
GG,-1.09,-6.09,-15.5
CC,-1.09,-6.09,-15.5
EA,0.49,20.73,62.6
UE,0.49,20.73,62.6
AE,0.48,20.20,61.0
EU,0.48,20.20,61.0
EC,0.63,21.07,63.2
GE,0.63,21.07,63.2
CE,0.43,18.09,54.7
EG,0.43,18.09,54.7
"""

# Table 3 — phosphorothioate / phosphodiester (PS/PO)
ps_po_table = """
AA,-0.52,-5.81,-16.4
AU,-0.42,-5.64,-16.2
AC,-0.88,-8.66,-24.1
AG,-0.71,-6.13,-16.8
UA,-0.30,-3.86,-11.0
UU,-0.49,-5.87,-16.7
UC,-0.64,-6.13,-17.0
UG,-0.77,-7.28,-20.2
CA,-0.82,-7.23,-19.8
CU,-0.66,-6.30,-17.5
CC,-0.91,-5.57,-14.4
CG,-1.21,-8.07,-21.2
GA,-0.85,-7.92,-21.9
GU,-0.61,-5.46,-15.0
GC,-1.15,-6.63,-17.0
GG,-1.09,-7.75,-20.6
EA,-0.57,13.64,44.0
AE,-0.54,15.07,48.3
EU,-0.59,14.87,47.8
UE,-0.56,14.72,47.3
EC,-0.58,10.31,33.7
CE,-0.55,10.49,34.2
EG,-0.35,16.13,51.0
GE,-0.44,14.65,46.7
"""


def _gibbs_37(table: str) -> "pd.Series":
    df = pd.read_csv(StringIO(table), header=None, names=["nucleotide", "G_50", "H", "S"])
    df["G_37"] = df["H"] - _BODY_TEMPERATURE_K * (df["S"] / 1000.0)
    return df.set_index("nucleotide")["G_37"]


po_po_g37 = _gibbs_37(po_po_table)
ps_po_g37 = _gibbs_37(ps_po_table)

# PS backbone contribution per dinucleotide at 37 C (kcal/mol).
ps_delta_g37 = (ps_po_g37 - po_po_g37).to_dict()
