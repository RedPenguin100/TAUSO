import re
import numpy as np
import primer3
import codonbias as cb
import pandas as pd

from numba import njit
from scipy.stats import entropy
from collections import Counter

from ..util import get_antisense

from ..algorithms.suffix_array import longest_prefix

from primer3 import calc_hairpin
from collections import defaultdict

def hairpin_dG_energy(seq: str):
    """
    Returns the raw ΔG (Gibbs free energy) of predicted hairpin structure,
    divided by sequence length.

    This value is NOT normalized to [0,1].
    Positive values suggest unstable/no structure.
    Negative values indicate stronger/stable hairpins that may interfere with ASO activity.
    """
    hairpin = calc_hairpin(seq)

    if not hairpin.structure_found:
        return 0
    return hairpin.dg if len(seq) > 0 else 0


@njit
def _is_palindrome(seq: str) -> bool:
    final = ''
    test = seq[::-1]
    for i in range(len(test)):
        if test[i] == 'A':
            final += 'T'
        elif test[i] == 'C':
            final += 'G'
        elif test[i] == 'G':
            final += 'C'
        elif test[i] == 'T':
            final += 'A'
    return final == seq


@njit
def palindromic_fraction(seq: str, l: int) -> float:
    count = 0
    for n in range(len(seq) - l + 1):
        curr_seq = seq[n:n + l]
        count += _is_palindrome(curr_seq)
    return count / len(seq)


@njit
def homooligo_count(seq: str) -> float:
    seq += '$'
    tot_count = 0
    curr_seq = ''
    n = 0
    while n in range(len(seq) - 1):
        while seq[n] == seq[n + 1]:
            curr_seq += seq[n]
            n += 1
        if len(curr_seq) > 1:
            tot_count += len(curr_seq) + 1
        n += 1
        curr_seq = ''
    return tot_count / len(seq)


def compute_ENC(seq: str) -> float:
    """
    Returns normalized ENC in [0,1], or NaN if the input sequence is empty/invalid.
    """
    if not isinstance(seq, str) or seq.strip() == "":
        return np.nan

    enc = cb.scores.EffectiveNumberOfCodons(bg_correction=True)
    enc_score = enc.get_score(seq)

    if enc_score is None or not np.isfinite(enc_score):
        return np.nan

    # Normalize ENC from [20..61] to [0..1]
    return (enc_score - 20.0) / (61.0 - 20.0)


def seq_entropy(seq: str) -> float:
    freqs = [seq.count(base) / len(seq) for base in "ACGT"]
    return entropy(freqs) / 2


def count_g_runs(seq: str, min_run_length: int = 4) -> float:
    """
    Calculates the fraction of the sequence that contains G-runs
    of at least 'min_run_length' (like 'GGGG').
    Normalized by the sequence length.
    """
    if len(seq) == 0:
        return 0.0
    seq = seq.upper()
    count = 0
    i = 0
    while i < len(seq):
        if seq[i] == 'G':
            run_length = 1
            while i + 1 < len(seq) and seq[i + 1] == 'G':
                run_length += 1
                i += 1
            if run_length >= min_run_length:
                count += 1
        i += 1
    return count / len(seq)


def hairpin_score(seq: str, min_overlap: int = 4) -> float:
    """
    Estimates the potential of the sequence to form a hairpin structure
    by checking how many small subsequences appear in its reverse complement.
    """
    seq = seq.upper()
    antisense = get_antisense(seq)
    matches = 0
    for i in range(len(seq) - min_overlap + 1):
        sub = seq[i:i + min_overlap]
        if sub in antisense:
            matches += 1
    return matches / len(seq)

def gc_skew(seq: str) -> float:
    """
    Computes GC skew = (G - C) / (G + C)
    A measure of strand asymmetry in G/C content, can affect hybridization and folding.
    """
    seq = seq.upper()
    G_counts = seq.count("G")
    C_counts = seq.count("C")
    if G_counts + C_counts == 0:
        return  0.0
    return (G_counts - C_counts) / (G_counts + C_counts)

@njit
def get_gc_content(seq: str) -> float:
    gc_count = 0
    for i in range(len(seq)):
        if seq[i] in "GCgc":
            gc_count +=1

    return gc_count / len(seq)

def gc_content_3prime_end(aso_sequence: str, window: int = 5) -> float:
    """Calculate the GC content at the 3' end of the ASO sequence."""
    if len(aso_sequence) < window:
        return 0.0
    three_prime_end = aso_sequence[-window:]
    gc_count = three_prime_end.count('G') + three_prime_end.count('C')
    return gc_count / window

@njit
def gc_skew_ends(seq: str, window: int = 5) -> float:
    """
    Calculates the GC-content difference between the 5' and 3' ends of the sequence.
    Measures thermodynamic asymmetry between ends.
    """
    seq = seq.upper()
    if len(seq) < 2 * window:
        return 0.0
    start = seq[:window]
    end = seq[-window:]
    gc_5 = start.count("G") + start.count("C")
    gc_3 = end.count("G") + end.count("C")
    return (gc_5 - gc_3) / window

def dispersed_repeats_score(seq, min_unit=2, max_unit=6):
    """
    Counts motifs (2–6 nt) that appear more than once, even if not consecutive.
    Helps detect internal similarity and potential self-binding regions.
    """
    unit_counter = Counter()
    for unit_len in range(min_unit, max_unit + 1):
        for i in range(len(seq) - unit_len + 1):
            unit = seq[i:i+unit_len]
            unit_counter[unit] += 1
    score = sum(count - 1 for count in unit_counter.values() if count > 1)
    return score / len(seq)



@njit
def at_skew(seq: str) -> float:
    """
    Calculates AT skew = (A - T) / (A + T)
    A measure of asymmetry in A/T content, can affect flexibility and binding dynamics.
    """
    seq = seq.upper()
    A_counts = seq.count("A")
    T_counts = seq.count("T")
    if A_counts + T_counts == 0:
        return 0.0  # avoid division by zero
    return (A_counts - T_counts) / (A_counts + T_counts)


def toxic_motif_count(aso_sequence, motifs=['UGU', 'GGTGG', 'TGGT', 'GGGU']) -> float:
    """
    Counts the number of toxic motif appearances in the ASO.
    Returns normalized count (0–1) based on max possible motif hits.

    Parameters:
        aso_sequence (str): DNA or RNA ASO sequence
        motifs (list of str): known toxic motifs (excluding 'GGGG' to avoid overlap with G-runs)

    Returns:
        float: normalized toxic motif count
    """
    sequence = aso_sequence.upper()
    total = 0
    for motif in motifs:
        total += len(re.findall(motif, sequence))

    max_possible = len(sequence)  # conservative upper bound
    return min(total / max_possible, 1.0)


def nucleotide_diversity(seq: str) -> float:
    # checking the nucleotide diversity of the ASO sequence and normalize it by the
    # max value 16
    nucs = [seq[i:i+2] for i in range(len(seq)-1)]
    unique = set(nucs)
    return len(unique) / 16

def stop_codon_count(seq: str, codons=('TAA', 'TAG', 'TGA')) -> float:
    """
    Counts occurrences of stop codons (TAA, TAG, TGA) in the sequence.
    Returns a normalized value (count per length).
    Returns:
        float: normalized count of stop codons
    """
    seq = seq.upper()
    count = sum(seq.count(codon) for codon in codons)
    return count / len(seq)

def tandem_repeats_score(seq: str, min_unit=2, max_unit=6) -> float:
    """
    Calculates how many short motifs (2–6 nt) repeat consecutively (tandem repeats).
    Useful for detecting repetitive structures that may form hairpins or reduce specificity.
    """
    score = 0
    for unit_len in range(min_unit, max_unit + 1):
        for i in range(len(seq) - unit_len * 2 + 1):
            unit = seq[i:i+unit_len]
            repeat_count = 1
            j = i + unit_len
            while j + unit_len <= len(seq) and seq[j:j+unit_len] == unit:
                repeat_count += 1
                j += unit_len
            if repeat_count >= 2:
                score += repeat_count - 1
    return score/len(seq)

def flexible_dinucleotide_fraction(seq: str) -> float:
    """
    Computes the fraction of flexible dinucleotides (AT and TA) in the sequence.
    These dinucleotides are less stable and can indicate structurally flexible regions.
    Args:
        seq (str): DNA sequence (assumed to be uppercase A/C/G/T)
    Returns:
        float: Fraction of AT or TA dinucleotides, normalized by sequence length
    """
    if len(seq) < 2:
        return 0.0
    seq = seq.upper()
    count = 0
    for i in range(len(seq)-1):
        pair = seq[i:i+2]
        if pair in ['AT','TA']:
            count += 1
    return count/(len(seq)-1)

def hairpin_tm(seq: str) -> float:
    """
    Feature: hairpin_tm
    Melting temperature (Tm) of the predicted hairpin structure.
    Higher Tm means the structure is more stable at physiological temperatures.

    Returns:
        float: Tm of the hairpin
    """
    if len(seq) == 0:
        return 0.
    hairpin = primer3.calc_hairpin(seq)
    if not hairpin.structure_found:
        return 0.
    return hairpin.tm


####################################################################
def calculate_chimera_ars(suffix_array, target_sequence, step_size):
    longest_prefix_lengths = []

    for start_index in range(1, len(target_sequence), step_size):
        prefix, _ = longest_prefix(target_sequence[start_index:], suffix_array)
        longest_prefix_lengths.append(len(prefix))

    chimera_ars_score = np.mean(longest_prefix_lengths)
    return chimera_ars_score
###################################################################################
def add_interaction_features(df : pd.DataFrame, feature_pairs: list[tuple[str, str]]) -> pd.DataFrame:
    """
    Given a DataFrame and a list of columns pairs (colA ,colB), create new cloumns named "colA*colB"
    containing element-wise product.
    df - DataFrame with your base features already computed
    feature_pairs : list of(str,str)
    """
    for col_a , col_b in feature_pairs:
        new_col = f"{col_a}*{col_b}"
        if new_col in df.columns:
            continue
        df[new_col] = df[col_a] * df[col_b]
    return df
#################################################################

def cg_dinucleotide_fraction(seq: str) -> float:
    """
    Calculate the fraction of 'CG' dinucleotides within a DNA sequence.

    A higher CG dinucleotide fraction suggests increased potential for immunogenicity
    issues when used as antisense oligonucleotides (ASOs).

    Args:
        seq (str): DNA sequence (A, C, G, T)

    Returns:
        float: CG dinucleotide fraction (normalized between 0 and 1)
    """
    seq = seq.upper()
    cg_count = 0
    for i in range (len(seq)-1):
        dinucleotide = seq[i:i+2]
        if dinucleotide == 'CG':
            cg_count += 1
    total_possible_pairs = len(seq) - 1
    if total_possible_pairs == 0:
        return 0.0
    cg_fraction = cg_count / total_possible_pairs
    return cg_fraction
########################################################################

def poly_pyrimidine_stretch(seq: str, min_run_length: int = 4) -> float:
    """
    Calculates the fraction of sequence containing poly-pyrimidine stretches (C and T bases).

    Long consecutive pyrimidine runs may cause undesired secondary structures
    or non-specific binding of antisense oligonucleotides (ASOs).

    Args:
        seq (str): DNA sequence (A/C/G/T)
        min_run_length (int): Minimum length of pyrimidine run considered problematic (default=4)

    Returns:
        float: Normalized poly-pyrimidine stretch score between 0 and 1
    """
    seq = seq.upper()
    pyrimidines = "CT"
    stretch_count = 0
    i = 0

    while i < len(seq):
        run_length = 0
        while i < len(seq) and seq[i] in pyrimidines:
            run_length += 1
            i += 1

        if run_length >= min_run_length:
            stretch_count += 1

        if run_length == 0:
            i += 1

    return stretch_count / len(seq) if len(seq) > 0 else 0.0

############################################################################

def dinucleotide_entropy(seq: str) -> float:
    """
    Calculates normalized Shannon entropy (0–1) based on dinucleotide frequencies in the sequence.

    High entropy = greater structural complexity and diversity.
    Low entropy = repetitive or predictable dinucleotide patterns.

    Args:
        seq (str): DNA sequence

    Returns:
        float: Normalized entropy (0 to 1)
    """
    seq = seq.upper()
    if len(seq) < 2:
        return 0.0  # too short to form any dinucleotide

    dinucleotides = [seq[i:i+2] for i in range(len(seq)-1)]
    freq = pd.Series(dinucleotides).value_counts(normalize=True)
    raw_entropy = entropy(freq, base=2)

    return raw_entropy / 4  # normalization to range [0, 1]
##############################################################################
def gc_block_length(seq):
    """
    find the length of the longest cosecutive G/C block in the sequence
    """
    seq = seq.upper()
    max_len = 0
    curr_len = 0
    for base in seq:
        if base in "GC":
            curr_len += 1
            max_len = max(max_len , curr_len)
        else:
            curr_len = 0
    return max_len

################################################################################
def purine_content(seq):
    """
    Calculates the fraction of purine bases (A and G) in the sequence.
    Purine-rich sequences may be more stable and bind better to RNA targets.
    """
    seq = seq.upper()
    count = seq.count("A") + seq.count("G")
    if len(seq) == 0:
        return 0.0
    return count / len(seq)
#################################################################   
def Niv_ENC(seq: str, strict: bool = False) -> float:
    """
    Calculates the Effective Number of Codons (ENC) for a DNA sequence.
    
    If strict=True, uses the original Wright (1990) formula strictly,
    requiring all four F-values (F2, F3, F4, F6). Otherwise, uses only the 
    available families to compute a partial ENC approximation.
    
    Args:
        seq (str): DNA sequence (assumed uppercase A/C/G/T)
        strict (bool): Whether to enforce full Wright formula (default: False)
        
    Returns:
        float: Normalized ENC in [0, 1] (0 = max bias, 1 = no bias)
    """
    seq = seq.upper()
    seq = seq[:len(seq) - (len(seq) % 3)]  # Trim to full codons
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    codon_counts = defaultdict(int)
    for codon in codons:
        codon_counts[codon] += 1

    FAMILY_GROUPS = {
        2: [['TTT', 'TTC'], ['TAT', 'TAC'], ['CAT', 'CAC'], ['CAA', 'CAG'],
            ['AAT', 'AAC'], ['AAA', 'AAG'], ['GAT', 'GAC'], ['GAA', 'GAG'], ['TGT', 'TGC']],
        3: [['ATT', 'ATC', 'ATA']],
        4: [['CCT', 'CCC', 'CCA', 'CCG'], ['ACT', 'ACC', 'ACA', 'ACG'],
            ['GCT', 'GCC', 'GCA', 'GCG'], ['GTT', 'GTC', 'GTA', 'GTG'], ['GGT', 'GGC', 'GGA', 'GGG']],
        6: [['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
            ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
            ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG']]
    }

    F_values = {}
    for k, families in FAMILY_GROUPS.items():
        F_list = []
        for family in families:
            counts = [codon_counts[c] for c in family]
            total = sum(counts)
            if total == 0:
                continue
            freqs = [count / total for count in counts]
            F_aa = sum(f ** 2 for f in freqs)
            F_list.append(F_aa)
        if F_list:
            F_values[k] = np.mean(F_list)

    try:
        weights = {2: 9, 3: 1, 4: 5, 6: 3}

        if strict:
            # Require all four F-values to compute full ENC
            if not all(k in F_values for k in [2, 3, 4, 6]):
                return 0.0  # Not enough info to compute strict ENC
            ENC = 2 + 9 / F_values[2] + 1 / F_values[3] + 5 / F_values[4] + 3 / F_values[6]
        else:
            # Use only available F-values (partial ENC)
            ENC = 2 + sum(weights[k] / F_values[k] for k in F_values)

        normalized_enc = (ENC - 20) / (61 - 20)
        return max(0.0, min(1.0, normalized_enc))

    except ZeroDivisionError:
        return 1.0  # fallback if unexpected division by 0

########################################################################################################
def at_rich_region_score(seq: str, min_run_length: int = 4) -> float:
    """
    Calculates the fraction of the sequence containing AT-rich stretches (A and T bases).

    Long consecutive A/T runs may reduce structural stability and impact ASO performance.

    Args:
        seq (str): DNA sequence (A/C/G/T)
        min_run_length (int): Minimum length of A/T run considered problematic (default = 4)

    Returns:
        float: Normalized AT-rich region score (0 to 1)
    """
    seq = seq.upper()
    at_bases = "AT"
    stretch_count = 0
    i = 0
    while i < len(seq):
        run_length = 0
        while i < len(seq) and seq[i] in at_bases:
            run_length += 1
            i += 1
        if run_length >= min_run_length:
            stretch_count += 1
        if run_length == 0:
            i += 1
    return stretch_count / len(seq) if len(seq) > 0 else 0.0
