"""TAUSO row -> HELM string, the notation OligoGym's featurizers read.

A TAUSO row carries the chemistry as three parallel strings: `aso_sequence` (one base per
residue), `chemical_pattern` (one sugar character per residue) and `ps_pattern` (one linkage
character per inter-residue bond, so one shorter). HELM writes the same thing as a dotted list
of `sugar(base)phosphate` monomers, the last residue carrying no trailing phosphate:

    RNA1{[cEt](G)[sp].[cEt](C)[sp].d(T)[sp].[cEt](A)}$$$$V2.0

The monomer vocabulary below is the one OligoGym's own datasets use, read off the HELM strings
in `oligogym/resources/pkg_dataset/`.
"""

# TAUSO chemical_pattern character -> HELM sugar monomer.
SUGAR_TO_HELM = {
    "M": "[moe]",   # 2'-O-methoxyethyl
    "C": "[cEt]",   # constrained ethyl
    "L": "[lna]",   # locked nucleic acid
    "d": "d",       # deoxyribose
    "o": "m",       # 2'-O-methyl
    "f": "[fl2r]",  # 2'-fluoro
}

# TAUSO ps_pattern character -> HELM phosphate monomer.
LINKAGE_TO_HELM = {
    "*": "[sp]",  # phosphorothioate
    "d": "p",     # phosphodiester
}

HELM_TEMPLATE = "RNA1{{{monomers}}}$$$$V2.0"


def row_to_helm(sequence, chemical_pattern, ps_pattern):
    """HELM for one oligo. Raises ValueError if the three strings do not line up."""
    n = len(sequence)
    if len(chemical_pattern) != n:
        raise ValueError(f"chemical_pattern has {len(chemical_pattern)} chars for a {n}-mer: {chemical_pattern!r}")
    if len(ps_pattern) != n - 1:
        raise ValueError(f"ps_pattern has {len(ps_pattern)} chars for a {n}-mer: {ps_pattern!r}")

    monomers = []
    for i, (base, sugar_char) in enumerate(zip(sequence, chemical_pattern)):
        if sugar_char not in SUGAR_TO_HELM:
            raise ValueError(f"unknown sugar character {sugar_char!r} in {chemical_pattern!r}")
        if base not in "ACGTU":
            raise ValueError(f"unknown base {base!r} in {sequence!r}")
        phosphate = ""
        if i < n - 1:
            linkage_char = ps_pattern[i]
            if linkage_char not in LINKAGE_TO_HELM:
                raise ValueError(f"unknown linkage character {linkage_char!r} in {ps_pattern!r}")
            phosphate = LINKAGE_TO_HELM[linkage_char]
        monomers.append(f"{SUGAR_TO_HELM[sugar_char]}({base}){phosphate}")

    return HELM_TEMPLATE.format(monomers=".".join(monomers))


def add_helm_column(df, column="helm"):
    """Add a HELM column built from the sequence and the two chemistry patterns."""
    from tauso.data.consts import ASO_SEQUENCE, CHEMICAL_PATTERN, PS_PATTERN

    df = df.copy()
    df[column] = [
        row_to_helm(seq, chem, ps)
        for seq, chem, ps in zip(df[ASO_SEQUENCE], df[CHEMICAL_PATTERN], df[PS_PATTERN])
    ]
    return df
