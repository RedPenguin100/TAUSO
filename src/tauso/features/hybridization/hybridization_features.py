import logging

from ...common.modifications import check_pattern_length, get_longest_dna_gap
from ...util import (
    BODY_TEMPERATURE_C,
    DNA_BASES,
    celsius_to_kelvin,
    dna_to_rna,
    get_nucleotide_watson_crick,
    rna_to_dna,
)
from ..hybridization.exp_weights import DNA_RNA_DG37_WEIGHTS, PS_DELTA_DG37_WEIGHTS
from ..hybridization.weights.dna import DNA_DNA_WEIGHTS
from ..hybridization.weights.lna import LNA_DNA_WEIGHTS

logger = logging.getLogger(__name__)

# Nearest-neighbour sums below are keyed by the slices the loop already holds -- the two
# sequence characters, plus the two chemical-pattern characters where the sugar matters -- so a
# stack costs one dict lookup. The tables are built from the same per-dinucleotide helpers the
# loops fall back to, so the two paths cannot disagree.


def _dna_dna_increment(b1: str, b2: str):
    """(dH, dS) for one DNA/DNA dinucleotide over its Watson-Crick complement.

    Returns None when the table has no such stack. Raises on a base outside the DNA alphabet.
    """
    key = f"{b1}{b2}/{get_nucleotide_watson_crick(b1)}{get_nucleotide_watson_crick(b2)}"
    entry = DNA_DNA_WEIGHTS.get(key)
    if entry is None:
        logger.warning("Unknown key in weights table: %s", key)
        return None
    return entry["dH"], entry["dS"]


_DNA_DNA_BY_DINUCLEOTIDE = {b1 + b2: _dna_dna_increment(b1, b2) for b1 in DNA_BASES for b2 in DNA_BASES}


def _third_gen_increment(b1: str, b2: str, m1: str, m2: str, params: dict, letter: str):
    """(dH, dS) for one high-affinity-sugar stack, or None when the stack contributes nothing.

    A stack contributes only when it touches a ``letter`` sugar; pure-DNA ('dd') and any other
    sugar are skipped. Raises on a base outside the DNA alphabet.
    """
    if m1 == "d" and m2 == "d":
        return None
    if m1 == letter and m2 == letter:
        top = f"+{b1}+{b2}"
    elif m1 == "d" and m2 == letter:
        top = f"{b1}+{b2}"
    elif m1 == letter and m2 == "d":
        top = f"+{b1}{b2}"
    else:
        return None

    key = f"{top}/{get_nucleotide_watson_crick(b1)}{get_nucleotide_watson_crick(b2)}"
    entry = params.get(key)
    if entry is None:
        logger.warning("Unknown key in weights table: %s", key)
        return None
    return entry["dH"], entry["dS"]


_THIRD_GEN_TABLES: dict = {}


def _third_gen_table(params: dict, letter: str) -> dict:
    """``"<b1><b2><m1><m2>" -> (dH, dS) | None`` for the sugars that can contribute (``letter``
    and 'd'). Any other sugar misses the table and takes the ``_third_gen_increment`` path."""
    cached = _THIRD_GEN_TABLES.get(letter)
    if cached is not None and cached[0] is params:
        return cached[1]
    sugars = (letter, "d")
    table = {
        b1 + b2 + m1 + m2: _third_gen_increment(b1, b2, m1, m2, params, letter)
        for b1 in DNA_BASES
        for b2 in DNA_BASES
        for m1 in sugars
        for m2 in sugars
    }
    _THIRD_GEN_TABLES[letter] = (params, table)
    return table


def get_dna_rna_dg(seq: str) -> float:
    """Unmodified DNA/RNA hybrid dG (kcal/mol), summed 5'->3' over overlapping dinucleotides."""
    seq = dna_to_rna(seq)
    total = 0.0
    for i in range(len(seq) - 1):
        total += DNA_RNA_DG37_WEIGHTS[seq[i : i + 2]]
    return total


def get_ps_delta_dg(seq: str, ps_pattern: str) -> float:
    """Phosphorothioate backbone contribution relative to the DNA/RNA hybrid (kcal/mol)."""
    seq = dna_to_rna(seq)
    if not isinstance(ps_pattern, str) or len(ps_pattern) != len(seq) - 1:
        raise ValueError(f"ps_pattern must be a length-{len(seq) - 1} string for a {len(seq)}-mer, got {ps_pattern!r}")
    total = 0.0
    for i in range(len(seq) - 1):
        if ps_pattern[i] != "*":
            continue
        total += PS_DELTA_DG37_WEIGHTS[seq[i : i + 2]]
    return total


def get_ps_dna_rna_dg(seq: str, ps_pattern: str) -> float:
    """PS-modified DNA/RNA hybrid dG (kcal/mol): the DNA/RNA baseline plus the PS delta."""
    return get_dna_rna_dg(seq) + get_ps_delta_dg(seq, ps_pattern)


def get_dna_rna_dg_region(seq: str, chemical_pattern: str, region: str) -> float:
    """DNA/RNA hybrid dG (kcal/mol) for one region of a gapmer.

    A gapmer (and, approximately, a mixmer) has a central deoxy block that recruits RNase H
    flanked by two high-affinity wings that drive target binding. This splits the DNA/RNA dG
    into the 5' wing, the central gap, and the 3' wing. The gap is taken as the longest
    contiguous deoxy run in the chemical pattern; each overlapping dinucleotide is attributed
    to the region of its 5' base, so the three regions partition the full DNA/RNA dG.

    The sequence and chemical pattern are 5'->3' (the dataset's aso_sequence_5_to_3), so
    'wing5' is the 5'-terminal wing. Returns 0.0 for an empty region, NaN when the oligo has no
    chemical pattern, and raises when it has one of the wrong length.
    """
    check_pattern_length(seq, chemical_pattern)
    if not isinstance(chemical_pattern, str):
        return float("nan")

    gap_start, gap_end, gap_len = get_longest_dna_gap(chemical_pattern)
    if gap_len == 0:
        gap_start, gap_end = len(seq), len(seq)

    seq = dna_to_rna(seq)
    # Each dinucleotide belongs to the region of its 5' base, so a region is the half-open
    # index range [start, stop).
    bounds = {"wing5": (0, gap_start), "gap": (gap_start, gap_end), "wing3": (gap_end, len(seq))}
    start, stop = bounds.get(region, (0, 0))

    total = 0.0
    for i in range(start, min(stop, len(seq) - 1)):
        total += DNA_RNA_DG37_WEIGHTS[seq[i : i + 2]]
    return total


def calculate_3rd_gen_diff(seq, fmt, params, temp_c=BODY_TEMPERATURE_C, letter="L", region=None):
    """High-affinity sugar (LNA/cEt) delta dG (kcal/mol) summed 5'->3'.

    ``params`` holds nearest-neighbour increments keyed by the modified-strand dinucleotide
    ('+' marks a modified sugar) over the Watson-Crick complement. Only dinucleotides that
    touch a modified sugar contribute; pure-DNA ('dd') stacks are skipped.

    These increments are LNA-DNA/DNA values defined on top of the DNA/DNA baseline
    (SantaLucia & Hicks 2004), i.e. over a DNA target. No complete LNA-DNA/RNA nearest-
    neighbour set has been published, so they are used as a proxy for the A-form RNA-bound
    wing and kept as a standalone feature rather than summed onto the DNA/RNA baseline.

    ``region`` ('wing5'/'wing3') restricts the sum to one wing, split around the longest deoxy gap.
    """
    if len(seq) != len(fmt):
        return None

    # Detect the deoxy gap on the original pattern (lowercase 'd') before upper-casing.
    if region is not None:
        gap_start, gap_end, gap_len = get_longest_dna_gap(fmt)
        if gap_len == 0:
            gap_start, gap_end = len(fmt), len(fmt)

    seq = seq.upper()
    temp_k = celsius_to_kelvin(temp_c)
    table = _third_gen_table(params, letter)

    if region == "wing5":
        start, stop = 0, gap_start
    elif region == "wing3":
        start, stop = gap_end, len(seq)
    else:
        start, stop = 0, len(seq)

    total_dH = 0.0
    total_dS = 0.0
    for i in range(start, min(stop, len(seq) - 1)):
        key = seq[i : i + 2] + fmt[i : i + 2]
        if key in table:
            increment = table[key]
        else:
            increment = _third_gen_increment(seq[i], seq[i + 1], fmt[i], fmt[i + 1], params, letter)
        if increment is not None:
            total_dH += increment[0]
            total_dS += increment[1]

    return total_dH - (temp_k * (total_dS / 1000.0))


def is_single_sugar_mod(chemical_pattern, sugar):
    """True if the ASO is a pure gapmer of one high-affinity ``sugar``: built only from that sugar
    and DNA ('d'), with ``sugar`` present at least once. Mixmers or any other modified sugar make it
    False, so the single-sugar affinity model applies only when this is True."""
    if not isinstance(chemical_pattern, str) or sugar not in chemical_pattern:
        return False
    return all(code in (sugar, "d") for code in chemical_pattern)


def get_cet_dna_rna_dg(antisense, chemical_pattern):
    """cEt affinity re-referenced to the RNA target (kcal/mol): the DNA/RNA baseline plus the cEt
    increment, on the same RNA-target footing as the 2'-MOE MD terms (the increment is LNA-DNA/DNA,
    used as a proxy; see calculate_cet_delta). NaN when the oligo is built from nucleotides other
    than cEt or DNA (mixmers and any other modified sugar), or on length mismatch.
    """
    if not is_single_sugar_mod(chemical_pattern, "C"):
        return float("nan")
    cet = calculate_cet_delta(antisense, chemical_pattern)
    return float("nan") if cet is None else get_dna_rna_dg(antisense) + cet


def get_lna_dna_rna_dg(antisense, chemical_pattern):
    """LNA affinity re-referenced to the RNA target (kcal/mol): the DNA/RNA baseline plus the LNA
    increment (see calculate_lna_delta). NaN when the oligo is built from nucleotides other than
    LNA or DNA (mixmers and any other modified sugar), or on length mismatch.
    """
    if not is_single_sugar_mod(chemical_pattern, "L"):
        return float("nan")
    lna = calculate_lna_delta(antisense, chemical_pattern)
    return float("nan") if lna is None else get_dna_rna_dg(antisense) + lna


def get_cet_wing_dg(antisense, chemical_pattern, region):
    """cEt wing affinity re-referenced to the RNA target (kcal/mol): the DNA/RNA baseline over the
    wing plus the cEt increment over the same wing ('wing5'/'wing3'), on the same RNA-target footing
    as hybr_cet_dna_rna_dg. NaN when the oligo is built from nucleotides other than cEt or DNA
    (mixmers and any other modified sugar)."""
    if not is_single_sugar_mod(chemical_pattern, "C"):
        return float("nan")
    v = calculate_3rd_gen_diff(antisense, chemical_pattern, LNA_DNA_WEIGHTS, letter="C", region=region)
    return float("nan") if v is None else get_dna_rna_dg_region(antisense, chemical_pattern, region) + v


def get_lna_wing_dg(antisense, chemical_pattern, region):
    """LNA wing affinity re-referenced to the RNA target (kcal/mol): the DNA/RNA baseline over the
    wing plus the LNA increment over the same wing ('wing5'/'wing3'), on the same RNA-target footing
    as hybr_lna_dna_rna_dg. NaN when the oligo is built from nucleotides other than LNA or DNA
    (mixmers and any other modified sugar)."""
    if not is_single_sugar_mod(chemical_pattern, "L"):
        return float("nan")
    v = calculate_3rd_gen_diff(antisense, chemical_pattern, LNA_DNA_WEIGHTS, letter="L", region=region)
    return float("nan") if v is None else get_dna_rna_dg_region(antisense, chemical_pattern, region) + v


def calculate_lna_delta(antisense, chemical_pattern):
    """LNA high-affinity-sugar delta dG (kcal/mol).

    Increments are LNA-DNA/DNA nearest-neighbour values (McTigue 2004; Owczarzy 2011), i.e.
    the LNA strand paired against a DNA target; see calculate_3rd_gen_diff for the reference
    state and why no LNA-DNA/RNA set is used.
    """
    return calculate_3rd_gen_diff(antisense, chemical_pattern, LNA_DNA_WEIGHTS, letter="L")


def calculate_cet_delta(antisense, chemical_pattern):
    """cEt high-affinity-sugar delta dG (kcal/mol).

    Heuristic (not a mistake): cEt has no published cEt-specific nearest-neighbour set, so it
    reuses the LNA parameters on structural-homology grounds, matching the article's treatment.
    Same LNA-DNA/DNA reference state as calculate_lna_delta.
    """
    return calculate_3rd_gen_diff(antisense, chemical_pattern, LNA_DNA_WEIGHTS, letter="C")


def calculate_dna(antisense, temp_c=BODY_TEMPERATURE_C):
    """Unmodified DNA/DNA duplex dG (kcal/mol) summed 5'->3' (SantaLucia & Hicks 2004)."""
    seq = rna_to_dna(antisense)
    temp_k = celsius_to_kelvin(temp_c)

    total_dH = 0.0
    total_dS = 0.0
    for i in range(len(seq) - 1):
        pair = seq[i : i + 2]
        if pair in _DNA_DNA_BY_DINUCLEOTIDE:
            increment = _DNA_DNA_BY_DINUCLEOTIDE[pair]
        else:
            increment = _dna_dna_increment(seq[i], seq[i + 1])
        if increment is not None:
            total_dH += increment[0]
            total_dS += increment[1]

    return total_dH - (temp_k * (total_dS / 1000.0))
