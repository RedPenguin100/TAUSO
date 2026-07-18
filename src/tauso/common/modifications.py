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


def is_gapmer(pattern):
    """A pattern is a "real" gapmer iff it has both flanks AND a deoxy gap.

    MMMdddMMM -> yes  (wings flank a DNA gap)
    MMMMMMMMM -> no   (fully modified: no gap)
    ddddddddd -> no   (all-DNA: no wings)
    dddddMMMM -> no   (gap reaches the 5' end: no 5' wing)
    """
    if not isinstance(pattern, str) or not pattern:
        return False
    start, end, gap_len = get_longest_dna_gap(pattern)
    return gap_len > 0 and start > 0 and end < len(pattern)


# --- IDT order-string notation --------------------------------------------------------------------
#
# Per-sugar IDT residue codes, keyed by base. 2'-MOE has distinct 5'-terminal / internal / 3'-terminal
# vendor codes; DNA and LNA use one code at any position. IDT's 2'-MOE cytidine (/i2MOErC/) is already
# the 5-methyl-C form, so 2'-MOE wing cytosines are methylated by construction; IDT's LNA cytidine
# (+C) is likewise 5-methyl by construction. A DNA cytosine renders as IDT's 5-methyl-dC by default
# (internal /iMe-dC/, 5'/3'-terminal /5Me-dC/ / /3Me-dC/); `no_methyl_cytosine` in modification_string
# opts out to plain dC. cEt (chemical-pattern code `C`) is not an IDT catalogue product, so it is
# rejected rather than rendered.
_IDT_DNA = {b: b for b in "ACGT"}
_IDT_DNA_5MEC_5PRIME = {**_IDT_DNA, "C": "/5Me-dC/"}
_IDT_DNA_5MEC_INTERNAL = {**_IDT_DNA, "C": "/iMe-dC/"}
_IDT_DNA_5MEC_3PRIME = {**_IDT_DNA, "C": "/3Me-dC/"}
_IDT_LNA = {b: f"+{b}" for b in "ACGT"}
_IDT_MOE_5PRIME = {b: f"/52MOEr{b}/" for b in "ACGT"}
_IDT_MOE_INTERNAL = {b: f"/i2MOEr{b}/" for b in "ACGT"}
_IDT_MOE_3PRIME = {b: f"/32MOEr{b}/" for b in "ACGT"}

# CHEMICAL_PATTERN letter -> (5'-terminal, internal, 3'-terminal) base->code tables.
_IDT_SUGAR_CODES = {
    "d": (_IDT_DNA, _IDT_DNA, _IDT_DNA),
    "L": (_IDT_LNA, _IDT_LNA, _IDT_LNA),
    "M": (_IDT_MOE_5PRIME, _IDT_MOE_INTERNAL, _IDT_MOE_3PRIME),
}
# DNA sugar tables when cytosines are 5-methylated (a 2'-MOE wing C is already 5-methyl via /i2MOErC/).
_IDT_DNA_5MEC_CODES = (_IDT_DNA_5MEC_5PRIME, _IDT_DNA_5MEC_INTERNAL, _IDT_DNA_5MEC_3PRIME)

_PS_MARK = "*"  # phosphorothioate linkage in ps_pattern; any other character = phosphodiester


def _is_5_methyl_c(modification_string) -> bool:
    """Cytosines are 5-methyl by default; 'no_methyl_cytosine' in modification_string opts out."""
    norm = str(modification_string or "").lower().replace(" ", "_").replace("-", "_")
    return "no_methyl_cytosine" not in norm


def to_idt_notation(
    sequence: str, chemical_pattern: str, ps_pattern: str = None, modification_string: str = None
) -> str:
    """Render an ASO as an IDT-style order string, read 5'->3'.

    Each residue is emitted in its IDT code by sugar and base -- sugar from ``chemical_pattern``
    (``M``=2'-MOE, ``L``=LNA, ``d``=DNA) -- and a ``*`` is inserted between residues for
    every phosphorothioate linkage in ``ps_pattern`` (any other linkage character is phosphodiester
    and emits nothing). ``ps_pattern`` has one character per inter-nucleotide linkage
    (``len == len(sequence) - 1``); if omitted, an all-phosphorothioate backbone is assumed.

    DNA cytosines render as 5-methyl-dC (``/iMe-dC/``) by default; pass
    ``modification_string="no_methyl_cytosine"`` for plain dC. 2'-MOE and LNA wing cytosines are
    already the 5-methyl form (``/i2MOErC/``, ``+C``).

    Only DNA / 2'-MOE / LNA sugars and A/C/G/T bases are supported; anything else raises ValueError.
    cEt (chemical-pattern code ``C``) is rejected because IDT does not sell cEt phosphoramidites, so a
    cEt-containing ASO cannot be turned into an IDT order string."""
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

    dna_codes = _IDT_DNA_5MEC_CODES if _is_5_methyl_c(modification_string) else (_IDT_DNA, _IDT_DNA, _IDT_DNA)
    sugar_codes = {**_IDT_SUGAR_CODES, "d": dna_codes}

    residues = []
    for i, (base, sugar) in enumerate(zip(seq, pattern)):
        if base not in _IDT_DNA:
            raise ValueError(f"unsupported base {base!r} at position {i}")
        if sugar == "C":
            raise ValueError(
                f"cEt sugar (code 'C') at position {i} cannot be ordered from IDT (not a catalogue "
                "modification); remove the cEt residues or order this oligo from a cEt supplier."
            )
        if sugar not in sugar_codes:
            raise ValueError(f"unsupported sugar {sugar!r} at position {i} (expected one of M/L/d)")
        five, internal, three = sugar_codes[sugar]
        table = five if i == 0 else three if i == n - 1 else internal
        residues.append(table[base])

    out = residues[0] if residues else ""
    for i in range(1, n):
        out += (_PS_MARK if ps_pattern[i - 1] == _PS_MARK else "") + residues[i]
    return out
