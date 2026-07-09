import math
from collections import Counter

import pandas as pd
import primer3
import ViennaRNA as RNA
from Bio.SeqUtils import gc_fraction
from primer3 import calc_hairpin

from ...util import get_antisense


def ry_transition_fraction(seq: str) -> float:
    """
    Calculates the fraction of steps that alternate between Purines (R) and Pyrimidines (Y).
    """
    if len(seq) < 2:
        return 0.0

    seq = seq.upper()
    purines = set("AG")
    pyrimidines = set("CTU")
    transitions = 0

    for i in range(len(seq) - 1):
        b1, b2 = seq[i], seq[i + 1]
        if (b1 in purines and b2 in pyrimidines) or (b1 in pyrimidines and b2 in purines):
            transitions += 1

    return transitions / (len(seq) - 1)


def keto_amino_skew(seq: str) -> float:
    """
    Calculates Keto-Amino skew: (K - M) / (K + M).
    K (Keto) = G, T, U. M (Amino) = A, C.
    """
    seq = seq.upper()
    keto_count = sum(seq.count(b) for b in "GTU")
    amino_count = sum(seq.count(b) for b in "AC")

    total = keto_count + amino_count
    if total == 0:
        return 0.0

    return (keto_count - amino_count) / total


def terminal_gc_fraction(seq: str) -> float:
    """Fraction of the four terminal bases (first two + last two) that are G or C."""
    if len(seq) < 4:
        return 0.0

    terminals = (seq[:2] + seq[-2:]).upper()
    return sum(1 for b in terminals if b in "GC") / 4.0


def self_energy(seq: str) -> float:
    return float(primer3.calc_homodimer(seq).dg)


def internal_fold_rna(seq: str) -> float:
    """Self-structure ΔG (kcal/mol) of the ASO folded with RNA (Turner 2004) parameters.

    A gapmer's 2'-modified wings are ribose-like, so RNA nearest-neighbour energetics can fit the
    self-structure. The ASO is read 5'->3' with T mapped to U.
    """
    return RNA.fold(seq.upper().replace("T", "U"))[1]


def hairpin_dG_energy(seq: str):
    """
    Returns the raw ΔG (Gibbs free energy, primer3 units) of the predicted hairpin structure.

    Not normalized. Returns 0 when no structure is found; negative values indicate
    stronger/more stable hairpins that may interfere with ASO activity.
    """
    hairpin = calc_hairpin(seq)

    if not hairpin.structure_found:
        return 0
    return hairpin.dg if len(seq) > 0 else 0


_complement_map = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")


def is_palindrome(subseq: str) -> bool:
    return subseq.translate(_complement_map) == subseq[::-1]


def palindromic_fraction(seq: str, l: int) -> float:
    seq_len = len(seq)
    if seq_len < l:
        return 0.0

    # A generator expression inside sum() is highly optimized in Python's C backend
    count = sum(1 for n in range(seq_len - l + 1) if is_palindrome(seq[n : n + l]))

    return count / seq_len


def homooligo_count(seq: str) -> float:
    seq_len = len(seq)

    # Matches original behavior for an empty string (returns 0.0)
    if seq_len == 0:
        return 0.0

    tot_count = 0
    current_run = 1

    for i in range(1, seq_len):
        if seq[i] == seq[i - 1]:
            current_run += 1
        else:
            if current_run > 2:
                tot_count += current_run
            current_run = 1

    # Check the final run at the end of the sequence
    if current_run > 2:
        tot_count += current_run

    return tot_count / seq_len


def seq_entropy(seq: str) -> float:
    seq_len = len(seq)
    if seq_len == 0:
        return 0.0

    ent = 0.0
    # Native string count is fast, and we avoid creating intermediate lists
    for base in "ACGT":
        count = seq.count(base)
        if count > 0:
            p = count / seq_len
            ent -= p * math.log2(p)

    return ent / 2.0  # normalize to [0, 1]: max entropy over 4 bases is log2(4) = 2


def hairpin_score(seq: str, min_overlap: int = 4) -> float:
    """
    Estimates hairpin potential in O(N) time using set lookups.
    """
    n = len(seq)
    if n < min_overlap:
        return 0.0

    seq = seq.upper()

    # 1. Fast reverse complement using the helper function
    antisense = get_antisense(seq)

    # 2. Build a set of all min_overlap-mers in the antisense for O(1) lookup
    antisense_kmers = {antisense[i : i + min_overlap] for i in range(n - min_overlap + 1)}

    # 3. Count matches without scanning the whole string every time
    matches = sum(1 for i in range(n - min_overlap + 1) if seq[i : i + min_overlap] in antisense_kmers)

    return matches / n


def gc_skew(seq: str) -> float:
    """GC skew = (G - C) / (G + C): strand asymmetry in G/C content."""
    seq = seq.upper()
    G_counts = seq.count("G")
    C_counts = seq.count("C")
    if G_counts + C_counts == 0:
        return 0.0
    return (G_counts - C_counts) / (G_counts + C_counts)


def get_gc_content(seq: str) -> float:
    """GC fraction in [0, 1] via Biopython; 0.0 for an empty sequence."""
    if not seq:
        return 0.0
    return gc_fraction(seq)


def gc_content_3prime_end(aso_sequence: str, window: int = 5) -> float:
    """Calculate the GC content at the 3' end of the ASO sequence."""
    if len(aso_sequence) < window:
        return 0.0
    return get_gc_content(aso_sequence[-window:])


def gc_skew_ends(seq: str, window: int = 5) -> float:
    """GC-content difference between the 5' and 3' terminal windows of the sequence."""
    if len(seq) < 2 * window:
        return 0.0

    start = seq[:window]
    end = seq[-window:]

    # Count both upper and lower case directly on the tiny slices
    gc_5 = start.count("G") + start.count("C") + start.count("g") + start.count("c")
    gc_3 = end.count("G") + end.count("C") + end.count("g") + end.count("c")

    return (gc_5 - gc_3) / window


def dispersed_repeats_score(seq: str, min_unit: int = 2, max_unit: int = 6) -> float:
    """
    Counts motifs (2–6 nt) that appear more than once, even if not consecutive.
    Helps detect internal similarity and potential self-binding regions.
    """
    n = len(seq)
    if n < min_unit:
        return 0.0

    # A single generator expression creates all k-mers.
    # We pass this directly to Counter(), which consumes it at C-speeds.
    kmers = (seq[i : i + unit_len] for unit_len in range(min_unit, max_unit + 1) for i in range(n - unit_len + 1))

    counts = Counter(kmers)

    # Calculate the score. Using a generator here avoids building an intermediate list.
    score = sum(count - 1 for count in counts.values() if count > 1)

    return score / n


def at_skew(seq: str) -> float:
    """AT skew = (A - T) / (A + T): asymmetry in A/T content."""
    # Count both cases to avoid creating a new string in memory with .upper()
    A_counts = seq.count("A") + seq.count("a")
    T_counts = seq.count("T") + seq.count("t")

    total_AT = A_counts + T_counts
    if total_AT == 0:
        return 0.0  # avoid division by zero

    return (A_counts - T_counts) / total_AT


def dinucleotide_diversity(seq: str) -> float:
    # fraction of the 16 possible dinucleotides present in the ASO sequence
    nucs = [seq[i : i + 2] for i in range(len(seq) - 1)]
    unique = set(nucs)
    return len(unique) / 16


def tandem_repeats_score(seq: str, min_unit=2, max_unit=6) -> float:
    """
    Counts short motifs (2–6 nt) that repeat consecutively (tandem repeats), normalized by length.
    """
    if len(seq) == 0:
        return 0.0
    score = 0
    for unit_len in range(min_unit, max_unit + 1):
        for i in range(len(seq) - unit_len * 2 + 1):
            unit = seq[i : i + unit_len]
            repeat_count = 1
            j = i + unit_len
            while j + unit_len <= len(seq) and seq[j : j + unit_len] == unit:
                repeat_count += 1
                j += unit_len
            if repeat_count >= 2:
                score += repeat_count - 1
    return score / len(seq)


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
    for i in range(len(seq) - 1):
        pair = seq[i : i + 2]
        if pair in ["AT", "TA"]:
            count += 1
    return count / (len(seq) - 1)


def hairpin_tm(seq: str) -> float:
    """
    Feature: hairpin_tm
    Melting temperature (Tm) of the predicted hairpin structure.
    Higher Tm means the structure is more stable at physiological temperatures.

    Returns:
        float: Tm of the hairpin
    """
    if len(seq) == 0:
        return 0.0
    hairpin = primer3.calc_hairpin(seq)
    if not hairpin.structure_found:
        return 0.0
    return hairpin.tm


def add_interaction_features(df: pd.DataFrame, feature_pairs: list[tuple[str, str]]) -> pd.DataFrame:
    """
    Given a DataFrame and a list of columns pairs (colA ,colB), create new cloumns named "colA*colB"
    containing element-wise product.
    df - DataFrame with your base features already computed
    feature_pairs : list of(str,str)
    """
    for col_a, col_b in feature_pairs:
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
    for i in range(len(seq) - 1):
        dinucleotide = seq[i : i + 2]
        if dinucleotide == "CG":
            cg_count += 1
    total_possible_pairs = len(seq) - 1
    if total_possible_pairs == 0:
        return 0.0
    cg_fraction = cg_count / total_possible_pairs
    return cg_fraction


def ta_dinucleotide_fraction(seq: str) -> float:
    """
    Fraction of 'TA'/'UA' dinucleotides within a sequence, normalized by length.

    Args:
        seq (str): Nucleic acid sequence (handles both T and U)

    Returns:
        float: TA/UA dinucleotide fraction (normalized between 0 and 1)
    """
    seq = seq.upper()
    total_possible_pairs = len(seq) - 1

    if total_possible_pairs <= 0:
        return 0.0

    # Count both to safely handle either DNA or RNA inputs
    target_count = seq.count("TA") + seq.count("UA")

    return target_count / total_possible_pairs


def poly_pyrimidine_run_count(seq: str, min_run_length: int = 4) -> float:
    """
    Counts poly-pyrimidine (C/T) runs of length >= min_run_length, normalized by sequence length.

    Long consecutive pyrimidine runs may cause undesired secondary structures or non-specific binding
    of antisense oligonucleotides (ASOs). This is the NUMBER of such runs (not their length); the
    length of the longest run is captured separately by ``poly_pyrimidine_longest``.

    Args:
        seq (str): DNA sequence (A/C/G/T)
        min_run_length (int): Minimum length of pyrimidine run considered problematic (default=4)

    Returns:
        float: number of qualifying runs divided by sequence length
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


def dinucleotide_entropy(seq: str) -> float:
    """
    Calculates normalized Shannon entropy (0–1) based on dinucleotide frequencies.
    Optimized for speed using standard library modules.
    """
    seq = seq.upper()
    n = len(seq)

    if n < 2:
        return 0.0

    # Generator expression avoids building a large list in memory
    counts = Counter(seq[i : i + 2] for i in range(n - 1))

    total_dinucleotides = n - 1
    raw_entropy = 0.0

    # Calculate entropy manually to bypass SciPy overhead
    for count in counts.values():
        probability = count / total_dinucleotides
        raw_entropy -= probability * math.log2(probability)

    return raw_entropy / 4.0  # Normalize to [0, 1] (since max entropy for 16 states is log2(16) = 4)


##############################################################################
def gc_longest(seq: str) -> float:
    """Longest consecutive G/C run, normalized by sequence length."""
    return _longest_run_fraction(seq, "GC")


################################################################################
def purine_content(seq):
    """Fraction of purine bases (A and G) in the sequence."""
    seq = seq.upper()
    count = seq.count("A") + seq.count("G")
    if len(seq) == 0:
        return 0.0
    return count / len(seq)


########################################################################################################
def at_rich_run_count(seq: str, min_run_length: int = 4) -> float:
    """
    Counts AT-rich (A/T) runs of length >= min_run_length, normalized by sequence length.

    Long consecutive A/T runs may reduce structural stability and impact ASO performance. This is the
    NUMBER of such runs (not their length); the longest run is captured by ``at_rich_longest``.

    Args:
        seq (str): DNA sequence (A/C/G/T)
        min_run_length (int): Minimum length of A/T run considered problematic (default = 4)

    Returns:
        float: number of qualifying runs divided by sequence length
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


def _longest_run_fraction(seq: str, bases: str) -> float:
    """Length of the longest consecutive run drawn from ``bases``, normalized by sequence length."""
    seq = seq.upper()
    n = len(seq)
    if n == 0:
        return 0.0
    longest = current = 0
    for base in seq:
        if base in bases:
            current += 1
            longest = max(longest, current)
        else:
            current = 0
    return longest / n


def poly_pyrimidine_longest(seq: str) -> float:
    """Longest consecutive pyrimidine (C/T) run, normalized by sequence length.

    Complements ``poly_pyrimidine_run_count`` (which counts how many runs of length >= 4 occur): this
    captures the length of the single worst tract rather than how many tracts there are.
    """
    return _longest_run_fraction(seq, "CT")


def at_rich_longest(seq: str) -> float:
    """Longest consecutive A/T run, normalized by sequence length.

    Complements ``at_rich_run_count`` (a count of A/T runs of length >= 4): this captures the
    length of the single longest A/T tract.
    """
    return _longest_run_fraction(seq, "AT")
