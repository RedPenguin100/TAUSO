import numpy as np

DNA_DNA_WEIGHTS = {
    # Source: SantaLucia & Hicks, Annu. Rev. Biophys. Biomol. Struct. 2004
    # Table 1: Unified Nearest-Neighbor Parameters for DNA/DNA
    # Identity Pairs (Sequence 5'->3' / Complement 3'<-5')
    "AA/TT": {"dH": -7.6, "dS": -21.3},
    "TT/AA": {"dH": -7.6, "dS": -21.3},  # Symmetry
    "AT/TA": {"dH": -7.2, "dS": -20.4},
    # TA/AT is distinct from AT/TA!
    "TA/AT": {"dH": -7.2, "dS": -21.3},
    "CA/GT": {"dH": -8.5, "dS": -22.7},
    "TG/AC": {"dH": -8.5, "dS": -22.7},  # Symmetry
    "GT/CA": {"dH": -8.4, "dS": -22.4},
    "AC/TG": {"dH": -8.4, "dS": -22.4},  # Symmetry
    "CT/GA": {"dH": -7.8, "dS": -21.0},
    "AG/TC": {"dH": -7.8, "dS": -21.0},  # Symmetry
    "GA/CT": {"dH": -8.2, "dS": -22.2},
    "TC/AG": {"dH": -8.2, "dS": -22.2},  # Symmetry
    "CG/GC": {"dH": -10.6, "dS": -27.2},
    # GC/CG is distinct!
    "GC/CG": {"dH": -9.8, "dS": -24.4},
    "GG/CC": {"dH": -8.0, "dS": -19.9},
    "CC/GG": {"dH": -8.0, "dS": -19.9},  # Symmetry
}

INIT_GC, INIT_AT = 0.98, 1.03
"""SantaLucia helix initiation, kcal/mol."""

DNA_HAIRPIN_LOOP_PENALTY_RAW = (
    # Source: SantaLucia & Hicks, Annu. Rev. Biophys. Biomol. Struct. 2004
    # Hairpin loop initiation dG37, kcal/mol, indexed by the number of unpaired bases.
    # Transcribed from memory and not yet checked against the paper.
    # Index 0-2 pad the table so it can be indexed by loop length; a loop needs three bases.
    0.0,
    0.0,
    0.0,
    3.5,
    3.5,
    3.3,
    4.0,
    4.2,
    4.3,
    4.5,
    4.6,
    4.7,
    4.8,
    4.9,
    4.9,
    5.0,
    5.1,
    5.1,
    5.2,
    5.2,
    5.3,
    5.3,
    5.4,
    5.4,
    5.4,
    5.5,
    5.5,
    5.5,
    5.6,
    5.6,
    5.6,
)

DNA_HAIRPIN_LOOP_PENALTY = np.array(DNA_HAIRPIN_LOOP_PENALTY_RAW, dtype=np.float64)
