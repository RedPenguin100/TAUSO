import logging

from ...common.modifications import get_longest_dna_gap
from ...util import BODY_TEMPERATURE_C, celsius_to_kelvin, get_nucleotide_watson_crick
from ..hybridization.exp_weights import DNA_RNA_DG37_WEIGHTS, PS_DELTA_DG37_WEIGHTS
from ..hybridization.weights.dna import DNA_DNA_WEIGHTS
from ..hybridization.weights.lna import LNA_DNA_WEIGHTS

logger = logging.getLogger(__name__)


def _to_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def get_dna_rna_dg(seq: str) -> float:
    """Unmodified DNA/RNA hybrid dG (kcal/mol), summed 5'->3' over overlapping dinucleotides."""
    seq = _to_rna(seq)
    total = 0.0
    for i in range(len(seq) - 1):
        L, R = seq[i], seq[i + 1]
        total += DNA_RNA_DG37_WEIGHTS[L + R]
    return total


def get_ps_delta_dg(seq: str, ps_pattern: str) -> float:
    """Phosphorothioate backbone contribution relative to the DNA/RNA hybrid (kcal/mol)."""
    seq = _to_rna(seq)
    if not isinstance(ps_pattern, str) or len(ps_pattern) != len(seq) - 1:
        raise ValueError(f"ps_pattern must be a length-{len(seq) - 1} string for a {len(seq)}-mer, got {ps_pattern!r}")
    total = 0.0
    for i in range(len(seq) - 1):
        if ps_pattern[i] != "*":
            continue
        L, R = seq[i], seq[i + 1]
        total += PS_DELTA_DG37_WEIGHTS[L + R]
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
    'wing5' is the 5'-terminal wing. Returns 0.0 for an empty region, NaN if the pattern and
    sequence lengths disagree.
    """
    if not isinstance(chemical_pattern, str) or len(chemical_pattern) != len(seq):
        return float("nan")

    gap_start, gap_end, gap_len = get_longest_dna_gap(chemical_pattern)
    if gap_len == 0:
        gap_start, gap_end = len(seq), len(seq)

    seq = _to_rna(seq)
    total = 0.0
    for i in range(len(seq) - 1):
        if i < gap_start:
            base_region = "wing5"
        elif i < gap_end:
            base_region = "gap"
        else:
            base_region = "wing3"
        if base_region == region:
            total += DNA_RNA_DG37_WEIGHTS[seq[i] + seq[i + 1]]
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

    total_dH = 0.0
    total_dS = 0.0
    for i in range(len(seq) - 1):
        if region == "wing5" and not i < gap_start:
            continue
        if region == "wing3" and not i >= gap_end:
            continue
        b1, b2 = seq[i], seq[i + 1]
        m1, m2 = fmt[i], fmt[i + 1]

        if m1 == "d" and m2 == "d":
            continue
        if m1 == letter and m2 == letter:
            top = f"+{b1}+{b2}"
        elif m1 == "d" and m2 == letter:
            top = f"{b1}+{b2}"
        elif m1 == letter and m2 == "d":
            top = f"+{b1}{b2}"
        else:
            continue

        c1, c2 = get_nucleotide_watson_crick(b1), get_nucleotide_watson_crick(b2)
        key = f"{top}/{c1}{c2}"
        if key in params:
            total_dH += params[key]["dH"]
            total_dS += params[key]["dS"]
        else:
            logger.warning("Unknown key in weights table: %s", key)

    return total_dH - (temp_k * (total_dS / 1000.0))


def _sole_high_affinity_sugar(chemical_pattern, sugar):
    """True iff the pattern is built only from this high-affinity ``sugar`` and DNA ('d'), so the
    single-sugar affinity model applies. Anything else in the string -- a different modified sugar
    (mixmer) or any other/unknown code -- makes it False. ``sugar`` must appear at least once.
    """
    return (
        isinstance(chemical_pattern, str) and sugar in chemical_pattern and not (set(chemical_pattern) - {sugar, "d"})
    )


def get_cet_dna_rna_dg(antisense, chemical_pattern):
    """cEt affinity re-referenced to the RNA target (kcal/mol): the DNA/RNA baseline plus the cEt
    increment, on the same RNA-target footing as the 2'-MOE MD terms (the increment is LNA-DNA/DNA,
    used as a proxy; see calculate_cet_delta). NaN unless the oligo is built only from cEt and DNA
    (so mixmers and any other modified sugar are excluded), or on length mismatch.
    """
    if not _sole_high_affinity_sugar(chemical_pattern, "C"):
        return float("nan")
    cet = calculate_cet_delta(antisense, chemical_pattern)
    return float("nan") if cet is None else get_dna_rna_dg(antisense) + cet


def get_lna_dna_rna_dg(antisense, chemical_pattern):
    """LNA affinity re-referenced to the RNA target (kcal/mol): the DNA/RNA baseline plus the LNA
    increment (see calculate_lna_delta). NaN unless the oligo is built only from LNA and DNA (so
    mixed-chemistry oligos are excluded), or on length mismatch.
    """
    if not _sole_high_affinity_sugar(chemical_pattern, "L"):
        return float("nan")
    lna = calculate_lna_delta(antisense, chemical_pattern)
    return float("nan") if lna is None else get_dna_rna_dg(antisense) + lna


def get_cet_wing_dg(antisense, chemical_pattern, region):
    """cEt wing affinity re-referenced to the RNA target (kcal/mol): the DNA/RNA baseline over the
    wing plus the cEt increment over the same wing ('wing5'/'wing3'), on the same RNA-target footing
    as hybr_cet_dna_rna_dg. NaN unless the oligo is built only from cEt and DNA (mixmers / any other
    modified sugar excluded)."""
    if not _sole_high_affinity_sugar(chemical_pattern, "C"):
        return float("nan")
    v = calculate_3rd_gen_diff(antisense, chemical_pattern, LNA_DNA_WEIGHTS, letter="C", region=region)
    return float("nan") if v is None else get_dna_rna_dg_region(antisense, chemical_pattern, region) + v


def get_lna_wing_dg(antisense, chemical_pattern, region):
    """LNA wing affinity re-referenced to the RNA target (kcal/mol): the DNA/RNA baseline over the
    wing plus the LNA increment over the same wing ('wing5'/'wing3'), on the same RNA-target footing
    as hybr_lna_dna_rna_dg. NaN unless the oligo is built only from LNA and DNA (mixed-chemistry
    oligos excluded)."""
    if not _sole_high_affinity_sugar(chemical_pattern, "L"):
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
    seq = antisense.upper().replace("U", "T")
    temp_k = celsius_to_kelvin(temp_c)

    total_dH = 0.0
    total_dS = 0.0
    for i in range(len(seq) - 1):
        b1, b2 = seq[i], seq[i + 1]
        c1, c2 = get_nucleotide_watson_crick(b1), get_nucleotide_watson_crick(b2)
        key = f"{b1}{b2}/{c1}{c2}"
        if key in DNA_DNA_WEIGHTS:
            total_dH += DNA_DNA_WEIGHTS[key]["dH"]
            total_dS += DNA_DNA_WEIGHTS[key]["dS"]
        else:
            logger.warning("Unknown key in weights table: %s", key)

    return total_dH - (temp_k * (total_dS / 1000.0))
