import logging

from ...common.modifications import get_longest_dna_gap
from ...util import BODY_TEMPERATURE_C, celsius_to_kelvin, get_nucleotide_watson_crick
from ..hybridization.exp_weights import ps_delta_g37
from ..hybridization.weights.dna import DNA_DNA_WEIGHTS
from ..hybridization.weights.lna import LNA_DNA_WEIGHTS

logger = logging.getLogger(__name__)

# DNA/RNA hybrid nearest-neighbour free energies, dG37 in kcal/mol.
# Source: Sugimoto et al., Biochemistry 1995. https://pubs.acs.org/doi/10.1021/bi00035a029
# Tabulated values are dG37 * 100 (the integers below); they are divided to kcal/mol once.
_DNA_RNA_DG37_CENTI = {
    "UU": -100,
    "GU": -210,
    "CU": -180,
    "AU": -90,
    "UG": -90,
    "GG": -210,
    "CG": -170,
    "AG": -90,
    "UC": -130,
    "GC": -270,
    "CC": -290,
    "AC": -110,
    "UA": -60,
    "GA": -150,
    "CA": -160,
    "AA": -20,
}
DNA_RNA_DG37 = {pair: centi / 100.0 for pair, centi in _DNA_RNA_DG37_CENTI.items()}


def _to_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def get_dna_rna_dg(seq: str) -> float:
    """Unmodified DNA/RNA hybrid dG (kcal/mol), summed 5'->3' over overlapping dinucleotides."""
    rna = _to_rna(seq)
    total = sum(DNA_RNA_DG37[rna[i : i + 2]] for i in range(len(rna) - 1))
    return round(total, 3)


def get_ps_delta_dg(seq: str) -> float:
    """Phosphorothioate backbone contribution (kcal/mol); positive = PS destabilises the duplex."""
    rna = _to_rna(seq)
    total = sum(ps_delta_g37[rna[i : i + 2]] for i in range(len(rna) - 1))
    return round(total, 3)


def get_ps_dna_rna_dg(seq: str) -> float:
    """PS-modified DNA/RNA hybrid dG (kcal/mol): the DNA/RNA baseline plus the PS backbone contribution."""
    return round(get_dna_rna_dg(seq) + get_ps_delta_dg(seq), 3)


def get_dna_rna_dg_region(seq: str, chemical_pattern: str, region: str) -> float:
    """DNA/RNA hybrid dG (kcal/mol) restricted to one region of a gapmer.

    The DNA gap is the longest run of 'd' in the chemical pattern; the 5' and 3'
    wings flank it. Each overlapping dinucleotide is attributed to the region of
    its 5' base, so the three regions partition the full DNA/RNA dG without
    double counting. Returns 0.0 when the requested region is empty.
    """
    if not isinstance(chemical_pattern, str) or len(chemical_pattern) != len(seq):
        return float("nan")

    gap_start, gap_end, gap_len = get_longest_dna_gap(chemical_pattern)
    if gap_len == 0:
        gap_start, gap_end = len(seq), len(seq)

    rna = _to_rna(seq)
    total = 0.0
    for i in range(len(rna) - 1):
        if i < gap_start:
            base_region = "wing5"
        elif i < gap_end:
            base_region = "gap"
        else:
            base_region = "wing3"
        if base_region == region:
            total += DNA_RNA_DG37[rna[i : i + 2]]
    return round(total, 3)


def calculate_3rd_gen_diff(seq, fmt, params, temp_c=BODY_TEMPERATURE_C, letter="L"):
    """High-affinity sugar (LNA/cEt) delta dG (kcal/mol) summed 5'->3'.

    ``params`` holds nearest-neighbour increments keyed by the modified-strand
    dinucleotide ('+' marks a modified sugar) over the Watson-Crick complement.
    Only dinucleotides that touch a modified sugar contribute; pure-DNA ('dd')
    stacks are skipped.
    """
    seq = seq.upper()
    fmt = fmt.upper()

    if len(seq) != len(fmt):
        return None

    temp_k = celsius_to_kelvin(temp_c)

    total_dH = 0.0
    total_dS = 0.0

    for i in range(len(seq) - 1):
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

    total_dG = total_dH - (temp_k * (total_dS / 1000.0))
    return round(total_dG, 3)


def calculate_lna(antisense, chemical_pattern):
    return calculate_3rd_gen_diff(antisense, chemical_pattern, LNA_DNA_WEIGHTS, letter="L")


def calculate_cet(antisense, chemical_pattern):
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

    total_dG = total_dH - (temp_k * (total_dS / 1000.0))
    return round(total_dG, 3)
