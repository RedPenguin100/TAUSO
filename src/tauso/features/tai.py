import math
import re
import numpy as np

def _clean_to_dna(seq: str) -> str:
    """
    Normalize sequence to DNA: remove whitespace, uppercase, and convert U->T.
    """
    return re.sub(r'\s+', '', str(seq)).upper().replace('U', 'T')


def tai_weights(category: str) -> dict:
    """
    Returns tAI weights for the requested organism.
    - "sc": Saccharomyces cerevisiae (as in your original code, normalized).
    - "hm": Human, returned as a DNA-only dictionary (T only). Any U/T duplicates
            are merged with a deterministic rule to avoid key overwrites.
    """
    if category == "sc":
        weight_dict = dict()
        # @formatter:off
        weight_dict["AAA"] = 7; weight_dict["AAC"] = 10; weight_dict["AAG"] = 16.24; weight_dict["AAT"] = 5.9
        weight_dict["ACA"] = 4.0011; weight_dict["ACC"] = 7.92; weight_dict["ACG"] = 2.28; weight_dict["ACT"] = 11
        weight_dict["AGA"] = 11; weight_dict["AGC"] = 4; weight_dict["AGG"] = 4.52; weight_dict["AGT"] = 2.36
        weight_dict["ATA"] = 2.0013; weight_dict["ATC"] = 9.36; weight_dict["ATG"] = 10.64; weight_dict["ATT"] = 18.9

        weight_dict["CAA"] = 9; weight_dict["CAC"] = 7; weight_dict["CAG"] = 3.88; weight_dict["CAT"] = 4.13
        weight_dict["CCA"] = 10.0002; weight_dict["CCC"] = 1.44; weight_dict["CCG"] = 3.2; weight_dict["CCT"] = 2
        weight_dict["CGA"] = 0.0006; weight_dict["CGC"] = 4.32; weight_dict["CGG"] = 1; weight_dict["CGT"] = 6
        weight_dict["CTA"] = 3; weight_dict["CTC"] = 1; weight_dict["CTG"] = 0.96; weight_dict["CTT"] = 0.59

        weight_dict["GAA"] = 14; weight_dict["GAC"] = 16; weight_dict["GAG"] = 6.48; weight_dict["GAT"] = 9.44
        weight_dict["GCA"] = 5.0011; weight_dict["GCC"] = 7.92; weight_dict["GCG"] = 1.6; weight_dict["GCT"] = 11
        weight_dict["GGA"] = 3; weight_dict["GGC"] = 16; weight_dict["GGG"] = 2.96; weight_dict["GGT"] = 9.44
        weight_dict["GTA"] = 2.0014; weight_dict["GTC"] = 10.08; weight_dict["GTG"] = 2.64; weight_dict["GTT"] = 14

        weight_dict["TAA"] = 0; weight_dict["TAC"] = 8; weight_dict["TAG"] = 0; weight_dict["TAT"] = 4.72
        weight_dict["TCA"] = 3.0011; weight_dict["TCC"] = 7.92; weight_dict["TCG"] = 1.96; weight_dict["TCT"] = 11
        weight_dict["TGA"] = 0; weight_dict["TGC"] = 4; weight_dict["TGG"] = 6; weight_dict["TGT"] = 2.36
        weight_dict["TTA"] = 7; weight_dict["TTC"] = 10; weight_dict["TTG"] = 12.24; weight_dict["TTT"] = 5.9
        # @formatter:on
        # normalize to [0,1] as in your original
        m = max(weight_dict.values())
        return {k: v / m for k, v in weight_dict.items()}

    if category == "hm":
        # Original mixed-T/U map (your values). We will convert to DNA-only and merge collisions.
        raw = {
            'GCC': 0.0, 'GCG': 0.22727272727272727, 'GCU': 1.0, 'GCT': 0.0, 'GCA': 0.45454545454545453,
            'GGG': 0.21875, 'GGT': 0.0, 'GGU': 0.5, 'GGC': 1.0, 'GGA': 0.1875,
            'CCT': 0.0, 'CCU': 0.2, 'CCC': 0.0, 'CCG': 0.5, 'CCA': 1.0,
            'ACA': 0.36363636363636365, 'ACG': 0.2727272727272727, 'ACC': 0.0, 'ACT': 0.0, 'ACU': 1.0,
            'GTT': 0.0, 'GTC': 0.0, 'GTU': 1.0, 'GTA': 0.14285714285714285, 'GTG': 0.21428571428571427,
            'AGU': 0.18181818181818182, 'TCA': 0.2727272727272727, 'AGC': 0.36363636363636365, 'TCU': 1.0,
            'TCC': 0.0, 'TCT': 0.0, 'TCG': 0.22727272727272727, 'AGT': 0.0, 'AGG': 0.5909090909090909,
            'CGT': 0.0, 'CGA': 0.0, 'CGG': 0.09090909090909091, 'CGU': 0.5454545454545454, 'CGC': 0.0,
            'AGA': 1.0, 'TTG': 1.0, 'CTG': 0.1111111111111111, 'CTT': 0.0, 'CTU': 0.037037037037037035,
            'CTA': 0.2222222222222222, 'CTC': 0.07407407407407407, 'TTA': 0.5185185185185185, 'TTU': 0.5,
            'TTT': 0.0, 'TTC': 1.0, 'AAU': 0.5, 'AAT': 0.0, 'AAC': 1.0, 'AAA': 0.4, 'AAG': 1.0,
            'GAT': 0.0, 'GAC': 1.0, 'GAU': 0.5, 'GAA': 1.0, 'GAG': 0.6428571428571429,
            'CAU': 0.5, 'CAC': 1.0, 'CAT': 0.0, 'CAG': 0.6111111111111112, 'CAA': 1.0,
            'ATC': 0.0, 'ATG': 1.0, 'ATT': 0.0, 'ATA': 0.15384615384615385, 'ATU': 1.0,
            'TAC': 1.0, 'TAT': 0.0, 'TAU': 0.5, 'TAG': 0.0, 'TAA': 0.0,
            'TGG': 1.0, 'TGA': 0.0, 'TGC': 1.0, 'TGT': 0.0, 'TGU': 0.5
        }
        dna = {}
        for k, v in raw.items():
            kk = k.replace('U', 'T')  # normalize key to DNA
            # Merge collisions deterministically; here we keep the maximum weight seen.
            dna[kk] = max(v, dna.get(kk, float('-inf')))
        return dna

    raise ValueError("category must be 'sc' or 'hm'")


def calc_tAI(seq: str,
             weight_dictionary: dict,
             *,
             allow_stop_as_zero: bool = False) -> float:
    """
    Robust tAI:
    - cleans sequence to DNA (U->T),
    - requires length to be divisible by 3 (else returns NaN),
    - averages ONLY over codons that were actually included,
    - unknown codon -> NaN,
    - zero-weight codon policy is configurable via allow_stop_as_zero.
    """
    s = _clean_to_dna(seq)
    if not s or (len(s) % 3 != 0):
        return float('nan')

    log_sum = 0.0
    counted = 0

    for i in range(0, len(s), 3):
        codon = s[i:i+3]
        if len(codon) != 3:
            continue
        w = weight_dictionary.get(codon, np.nan)
        if np.isnan(w):
            # unseen/unknown codon -> no valid tAI
            return float('nan')

        if w == 0.0:
            if allow_stop_as_zero:
                # Strict: any stop/zero codon kills tAI
                return 0.0
            else:
                # Lenient: skip zeros but don't count them in the denominator
                continue

        log_sum += math.log(w)
        counted += 1

    return math.exp(log_sum / counted) if counted > 0 else float('nan')
