from numba.typed import Dict, List
from numba import njit, types
import numpy as np


def get_longer_string(s1: str, s2: str) -> str:
    return s1 if len(s1) >= len(s2) else s2

@njit
def get_purine_content(seq: str) -> float:
    """
    Calculates the fraction of purine bases (A and G) in the sequence.
    Purine-rich sequences may be more stable and bind better to RNA targets.
    """
    if len(seq) == 0:
        return 0.

    purine_count = 0
    for i in range(len(seq)):
        if seq[i] in "AaGg":
            purine_count +=1

    return purine_count / len(seq)



@njit
def get_nucleotide_watson_crick(nucleotide: chr) -> chr:
    if nucleotide == 'A':
        return 'T'
    if nucleotide == 'G':
        return 'C'
    if nucleotide == 'C':
        return 'G'
    if nucleotide == 'U' or nucleotide == 'T':
        return 'A'
    raise ValueError(f"Unknown codon {nucleotide}")


def _norm_rna_to_dna(seq: str) -> str:
    """Normalize RNA to DNA alphabet (U->T), uppercase, strip whitespace."""
    return str(seq).upper().replace('U', 'T').replace(' ', '').replace('\t', '').replace('\n', '')


def _to_str_seq(x) -> str:
    """
    Coerce sequence-like (list/np.array/Series) or string to a clean uppercase DNA string.
    Converts U->T and strips whitespace. Ensures slicing returns a plain string (avoids pandas iterable assignment).
    """
    if isinstance(x, str):
        s = x
    else:
        try:
            s = ''.join(list(x))
        except Exception:
            s = str(x)
    return s.replace(' ', '').replace('\t', '').replace('\n', '').replace('U', 'T').upper()

@njit
def get_antisense(sense: str) -> str:
    antisense = ''
    for n in range(len(sense) - 1, -1, -1):
        antisense += get_nucleotide_watson_crick(sense[n])
    return antisense


@njit
def get_codon_to_aa():
    codon_to_aa = Dict.empty(key_type=types.string, value_type=types.string)
    codon_to_aa['TTT'] = 'F'
    codon_to_aa['TTC'] = 'F'

    codon_to_aa['TTA'] = 'L'
    codon_to_aa['TTG'] = 'L'

    codon_to_aa['TCT'] = 'S'
    codon_to_aa['TCC'] = 'S'
    codon_to_aa['TCA'] = 'S'
    codon_to_aa['TCG'] = 'S'

    codon_to_aa['TAT'] = 'Y'
    codon_to_aa['TAC'] = 'Y'

    codon_to_aa['TGT'] = 'C'
    codon_to_aa['TGC'] = 'C'

    codon_to_aa['TGG'] = 'W'

    codon_to_aa['CTT'] = 'L'
    codon_to_aa['CTC'] = 'L'
    codon_to_aa['CTA'] = 'L'
    codon_to_aa['CTG'] = 'L'

    codon_to_aa['CCT'] = 'P'
    codon_to_aa['CCC'] = 'P'
    codon_to_aa['CCA'] = 'P'
    codon_to_aa['CCG'] = 'P'

    codon_to_aa['CAT'] = 'H'
    codon_to_aa['CAC'] = 'H'

    codon_to_aa['CAA'] = 'Q'
    codon_to_aa['CAG'] = 'Q'

    codon_to_aa['CGT'] = 'R'
    codon_to_aa['CGC'] = 'R'
    codon_to_aa['CGA'] = 'R'
    codon_to_aa['CGG'] = 'R'

    codon_to_aa['ATT'] = 'I'
    codon_to_aa['ATC'] = 'I'
    codon_to_aa['ATA'] = 'I'

    codon_to_aa['ATG'] = 'M'

    codon_to_aa['ACT'] = 'T'
    codon_to_aa['ACC'] = 'T'
    codon_to_aa['ACA'] = 'T'
    codon_to_aa['ACG'] = 'T'

    codon_to_aa['AAT'] = 'N'
    codon_to_aa['AAC'] = 'N'

    codon_to_aa['AAA'] = 'K'
    codon_to_aa['AAG'] = 'K'

    codon_to_aa['AGT'] = 'S'
    codon_to_aa['AGC'] = 'S'

    codon_to_aa['AGA'] = 'R'
    codon_to_aa['AGG'] = 'R'

    codon_to_aa['GTT'] = 'V'
    codon_to_aa['GTC'] = 'V'
    codon_to_aa['GTA'] = 'V'
    codon_to_aa['GTG'] = 'V'

    codon_to_aa['GCT'] = 'A'
    codon_to_aa['GCC'] = 'A'
    codon_to_aa['GCA'] = 'A'
    codon_to_aa['GCG'] = 'A'

    codon_to_aa['GAT'] = 'D'
    codon_to_aa['GAC'] = 'D'

    codon_to_aa['GAA'] = 'E'
    codon_to_aa['GAG'] = 'E'

    codon_to_aa['GGT'] = 'G'
    codon_to_aa['GGC'] = 'G'
    codon_to_aa['GGA'] = 'G'
    codon_to_aa['GGG'] = 'G'

    return codon_to_aa


@njit
def get_all_aa_codon_friends(input_aa):
    codon_to_aa = get_codon_to_aa()
    codon_friends = []
    for codon, aa in codon_to_aa.items():
        if aa == input_aa:
            codon_friends.append(codon)

    return codon_friends


@njit
def get_trna_dict():
    S_Cerevisiae_trna_anticodon_copy = Dict.empty(key_type=types.string, value_type=types.int64)

    S_Cerevisiae_trna_anticodon_copy['AGC'] = 11
    S_Cerevisiae_trna_anticodon_copy['TGC'] = 5
    S_Cerevisiae_trna_anticodon_copy['GCC'] = 16
    S_Cerevisiae_trna_anticodon_copy['CCC'] = 2
    S_Cerevisiae_trna_anticodon_copy['TCC'] = 3
    S_Cerevisiae_trna_anticodon_copy['AGG'] = 2
    S_Cerevisiae_trna_anticodon_copy['TGG'] = 10
    S_Cerevisiae_trna_anticodon_copy['AGT'] = 11
    S_Cerevisiae_trna_anticodon_copy['CGT'] = 1
    S_Cerevisiae_trna_anticodon_copy['TGT'] = 4
    S_Cerevisiae_trna_anticodon_copy['AAC'] = 14
    S_Cerevisiae_trna_anticodon_copy['CAC'] = 2
    S_Cerevisiae_trna_anticodon_copy['TAC'] = 2
    S_Cerevisiae_trna_anticodon_copy['GAA'] = 10
    S_Cerevisiae_trna_anticodon_copy['GTT'] = 10
    S_Cerevisiae_trna_anticodon_copy['CTT'] = 14
    S_Cerevisiae_trna_anticodon_copy['TTT'] = 7
    S_Cerevisiae_trna_anticodon_copy['GTC'] = 16
    S_Cerevisiae_trna_anticodon_copy['CTC'] = 2
    S_Cerevisiae_trna_anticodon_copy['TTC'] = 14
    S_Cerevisiae_trna_anticodon_copy['GTG'] = 7
    S_Cerevisiae_trna_anticodon_copy['CTG'] = 1
    S_Cerevisiae_trna_anticodon_copy['TTG'] = 9
    S_Cerevisiae_trna_anticodon_copy['AGA'] = 11
    S_Cerevisiae_trna_anticodon_copy['CGA'] = 1
    S_Cerevisiae_trna_anticodon_copy['TGA'] = 3
    S_Cerevisiae_trna_anticodon_copy['GCT'] = 4
    S_Cerevisiae_trna_anticodon_copy['ACG'] = 6
    S_Cerevisiae_trna_anticodon_copy['CCG'] = 1
    S_Cerevisiae_trna_anticodon_copy['CCT'] = 1
    S_Cerevisiae_trna_anticodon_copy['TCT'] = 11
    S_Cerevisiae_trna_anticodon_copy['GAG'] = 1
    S_Cerevisiae_trna_anticodon_copy['TAG'] = 3
    S_Cerevisiae_trna_anticodon_copy['CAA'] = 10
    S_Cerevisiae_trna_anticodon_copy['TAA'] = 7
    S_Cerevisiae_trna_anticodon_copy['AAT'] = 13
    S_Cerevisiae_trna_anticodon_copy['TAT'] = 2
    S_Cerevisiae_trna_anticodon_copy['CAT'] = 5
    S_Cerevisiae_trna_anticodon_copy['GTA'] = 8
    S_Cerevisiae_trna_anticodon_copy['GCA'] = 4
    S_Cerevisiae_trna_anticodon_copy['CCA'] = 6
    return S_Cerevisiae_trna_anticodon_copy


@njit
def get_all_codons():
    return ['TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC', 'TCA', 'TCG', 'TAT', 'TAC', 'TGT', 'TGC', 'TGG', 'CTT', 'CTC',
            'CTA', 'CTG', 'CCT', 'CCC', 'CCA', 'CCG', 'CAT', 'CAC', 'CAA', 'CAG', 'CGT', 'CGC', 'CGA', 'CGG', 'ATT',
            'ATC', 'ATA', 'ATG', 'ACT', 'ACC', 'ACA', 'ACG', 'AAT', 'AAC', 'AAA', 'AAG', 'AGT', 'AGC', 'AGA', 'AGG',
            'GTT', 'GTC', 'GTA', 'GTG', 'GCT', 'GCC', 'GCA', 'GCG', 'GAT', 'GAC', 'GAA', 'GAG', 'GGT', 'GGC', 'GGA',
            'GGG']


@njit
def trna_copy_number(anti_codon: str) -> int:
    if anti_codon % 3 != 0:
        raise ValueError(f"codon {anti_codon} must be divisible by 3")

    return get_trna_dict().get(anti_codon, 0)


# For ease A and I are interchangeable, even though Wobble is only with inosine
# Numbers are good only for S. Cerevisiae
@njit
def nucleotide_wobble(codon_third: chr, anticodon_first: chr) -> float:
    if codon_third == 'G':
        if anticodon_first == 'U' or anticodon_first == 'T':
            return 0.68
        if anticodon_first == 'C':
            return 0.
        return 1.
    if codon_third == 'U' or codon_third == 'T':
        if anticodon_first in ['I', 'A']:
            return 0.
        if anticodon_first == 'G':
            return 0.41
        return 1.
    if codon_third == 'C':
        if anticodon_first in ['I', 'A']:
            return 0.28
        if anticodon_first == 'G':
            return 0.
        return 1.
    if codon_third == 'A':
        if anticodon_first in ['I', 'A']:
            return 0.9999
        if anticodon_first == 'U' or anticodon_first == 'T':
            return 0.
        return 1.
    raise ValueError(f"Unknown codon {codon_third}")


@njit
def get_nucleotide_watson_crick(nucleotide):
    if nucleotide == 'A':
        return 'T'
    if nucleotide == 'G':
        return 'C'
    if nucleotide == 'C':
        return 'G'
    if nucleotide == 'U' or nucleotide == 'T':
        return 'A'
    raise ValueError(f"Unknown nucleotide {nucleotide}")


@njit
def is_nucleotide_wobble(codon_third: chr, anticodon_first: chr) -> bool:
    return 1. - nucleotide_wobble(codon_third, anticodon_first) > 1e-6  # Lowest Wobble yields 1e-4 so this is safe


@njit
def is_nucleotide_watson_crick(nucleotide: chr, anti_nucleotide: chr) -> bool:
    if nucleotide == 'U':
        nucleotide = 'T'

    return nucleotide == get_nucleotide_watson_crick(anti_nucleotide)


@njit
def is_single_trna_translating(trna: str, codon: str) -> bool:
    if not is_nucleotide_watson_crick(codon[0], trna[2]):
        return False
    if not is_nucleotide_watson_crick(codon[1], trna[1]):
        return False

    codon_third = codon[2]
    anticodon_first = trna[0]
    return is_nucleotide_wobble(codon_third, anticodon_first)


@njit
def calculate_tai(seq: str) -> float:
    if len(seq) % 3 != 0:
        raise ValueError(f"Sequence length {len(seq)} must be divisible by 3 ")

    trna_dict = get_trna_dict()

    codon_to_weight_dict = Dict.empty(key_type=types.string, value_type=types.float64)
    codon_to_aa = get_codon_to_aa()

    all_codons = get_all_codons()

    for codon in all_codons:
        weight_per_trna = []
        codon_friends = get_all_aa_codon_friends(codon_to_aa[codon])

        for trna, copy_number in trna_dict.items():
            # Wobble can't happen between different amino acids
            if get_antisense(trna) not in codon_friends:
                continue

            if is_single_trna_translating(trna, codon):
                s_ij = nucleotide_wobble(codon[2], trna[0])
                weight_per_trna.append((1 - s_ij) * copy_number)
        codon_to_weight_dict[codon] = sum(weight_per_trna)

    # we would like to exclude the last (stop) codon from the calculation
    weights = []

    for i in range(0, len(seq) - 3, 3):
        codon = seq[i: i + 3]
        weights.append(codon_to_weight_dict[codon])

    np_weights = np.array(weights) / np.max(weights)  # normalize so number is between 0 and 1
    return np.exp(np.mean(np.log(np_weights)))
