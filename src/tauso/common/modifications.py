import re


def get_longest_dna_gap(chemical_pattern: str, marker: str = "d") -> tuple[int, int, int]:
    """
    Finds the longest consecutive stretch of DNA markers.
    Returns (start_index, end_index, length). Returns (-1, -1, 0) if none found.
    """
    matches = list(re.finditer(f"{re.escape(marker)}+", chemical_pattern))
    if not matches:
        return -1, -1, 0

    # Find the longest match
    longest = max(matches, key=lambda m: m.end() - m.start())
    return longest.start(), longest.end(), longest.end() - longest.start()


def phosphorothioate_fraction(ps_pattern: str) -> float:
    """Fraction of inter-nucleotide bonds that are phosphorothioate ('*') in PS_PATTERN."""
    s = str(ps_pattern)
    return s.count("*") / len(s) if s else 0.0


def deoxy_sugar_fraction(chemical_pattern: str) -> float:
    """Fraction of sugars that are deoxyribose / DNA ('d') in CHEMICAL_PATTERN."""
    s = str(chemical_pattern)
    return s.count("d") / len(s) if s else 0.0


# --- IDT order-string notation --------------------------------------------------------------------
#
# Per-sugar IDT residue codes, keyed by base. 2'-MOE has distinct 5'-terminal / internal / 3'-terminal
# vendor codes; DNA, LNA and cEt use one code at any position. These specific codes should be
# confirmed against IDT's current modification list before ordering -- in particular cEt is not a
# standard IDT catalogue modification (the token here follows IDT's slash convention as a placeholder),
# and the 2'-MOE codes should be verified.
_IDT_DNA = {b: b for b in "ACGT"}
_IDT_LNA = {b: f"+{b}" for b in "ACGT"}
_IDT_CET = {b: f"/icEt{b}/" for b in "ACGT"}
_IDT_MOE_5PRIME = {b: f"/52MOEr{b}/" for b in "ACGT"}
_IDT_MOE_INTERNAL = {b: f"/i2MOEr{b}/" for b in "ACGT"}
_IDT_MOE_3PRIME = {b: f"/32MOEr{b}/" for b in "ACGT"}

# CHEMICAL_PATTERN letter -> (5'-terminal, internal, 3'-terminal) base->code tables.
_IDT_SUGAR_CODES = {
    "d": (_IDT_DNA, _IDT_DNA, _IDT_DNA),
    "L": (_IDT_LNA, _IDT_LNA, _IDT_LNA),
    "C": (_IDT_CET, _IDT_CET, _IDT_CET),
    "M": (_IDT_MOE_5PRIME, _IDT_MOE_INTERNAL, _IDT_MOE_3PRIME),
}

_PS_MARK = "*"  # phosphorothioate linkage in ps_pattern; any other character = phosphodiester


def to_idt_notation(sequence: str, chemical_pattern: str, ps_pattern: str = None) -> str:
    """Render an ASO as an IDT-style order string, read 5'->3'.

    Each residue is emitted in its IDT code by sugar and base -- sugar from ``chemical_pattern``
    (``M``=2'-MOE, ``C``=cEt, ``L``=LNA, ``d``=DNA) -- and a ``*`` is inserted between residues for
    every phosphorothioate linkage in ``ps_pattern`` (any other linkage character is phosphodiester
    and emits nothing). ``ps_pattern`` has one character per inter-nucleotide linkage
    (``len == len(sequence) - 1``); if omitted, an all-phosphorothioate backbone is assumed.

    Only DNA / 2'-MOE / cEt / LNA sugars and A/C/G/T bases are supported; anything else raises
    ValueError. The emitted modification codes should be confirmed against IDT's current catalogue
    before ordering (see the code tables above)."""
    seq = str(sequence).upper()
    pattern = str(chemical_pattern)
    n = len(seq)
    if len(pattern) != n:
        raise ValueError(f"chemical_pattern length {len(pattern)} != sequence length {n}")
    if ps_pattern is None:
        ps_pattern = _PS_MARK * (n - 1)
    ps_pattern = str(ps_pattern)
    if n > 1 and len(ps_pattern) != n - 1:
        raise ValueError(f"ps_pattern length {len(ps_pattern)} != number of linkages {n - 1}")

    residues = []
    for i, (base, sugar) in enumerate(zip(seq, pattern)):
        if base not in _IDT_DNA:
            raise ValueError(f"unsupported base {base!r} at position {i}")
        if sugar not in _IDT_SUGAR_CODES:
            raise ValueError(f"unsupported sugar {sugar!r} at position {i} (expected one of M/C/L/d)")
        five, internal, three = _IDT_SUGAR_CODES[sugar]
        table = five if i == 0 else three if i == n - 1 else internal
        residues.append(table[base])

    out = residues[0] if residues else ""
    for i in range(1, n):
        out += (_PS_MARK if ps_pattern[i - 1] == _PS_MARK else "") + residues[i]
    return out
