import numpy as np
import re
from collections import Counter


def compute_mod_fraction(pattern):
    """
    Returns the fraction of modified residues (non-'d' characters)
    in a given chemical pattern string.
    """
    pattern = str(pattern)
    total = len(pattern)
    if total == 0:
        return 0.0
    modified = sum(1 for c in pattern if c != 'd')
    return modified / total

##########################################################################
def compute_mod_type_count(pattern):
    """
    Returns the normalized number of unique modification characters (excluding 'd'),
    divided by the total number of modified residues. Returns 0 if no modifications.
    """
    pattern = str(pattern)
    mod_chars = [c for c in pattern if c != 'd']
    if not mod_chars:
        return 0.0

    unique_mods = set(mod_chars)
    return len(unique_mods) / len(mod_chars)


########################################################################

def compute_mod_5prime_run(pattern):
    """
    Returns the length of the longest consecutive stretch of modified residues
    (non-'d') starting from the 5' end (left side of the pattern).
    """
    pattern = str(pattern)
    run = 0
    for c in pattern:
        if c != 'd':
            run += 1
        else:
            break

    length = len(pattern)
    if length == 0:
        return 0.0
        
    return run / length
###########################################################################
def compute_mod_3prime_run(pattern):
    """
    Returns the length of the longest consecutive stretch of modified residues
    (non-'d') starting from the 3' end (right side of the pattern).
    """
    pattern = str(pattern)[::-1]  # reverse pattern to simulate scanning from 3'
    run = 0
    for c in pattern:
        if c != 'd':
            run += 1
        else:
            break

    length = len(pattern)
    if length == 0:
        return 0.0
        
    return run / length
###########################################################################
def compute_mod_min_distance_to_5prime(pattern):
    """
    Returns the normalized distance (0 to 1) of the first modified residue (non-'d')
    from the 5' end of the pattern. Returns -1.0 if no modifications are found or if pattern is empty.
    """
    pattern = str(pattern)
    n = len(pattern)
    if n == 0:
        return -1.0

    for i, c in enumerate(pattern):
        if c != 'd':
            return i / n

    return -1.0

###########################################################################
def compute_mod_min_distance_to_3prime(pattern):
    """
    Returns the normalized distance (0 to 1) of the first modified residue (non-'d')
    from the 3' end of the pattern. Returns -1.0 if no modifications are found.
    """
    pattern = str(pattern)
    n = len(pattern)
    if n == 0:
        return -1.0

    for i, c in enumerate(reversed(pattern)):
        if c != 'd':
            return i / n

    return -1.0

###########################################################################
def compute_mod_pos_std(pattern):
    """
    Returns the standard deviation of modified residue positions (non-'d') in the pattern.
    Returns -1 if there are no modifications.
    """


    pattern = str(pattern)
    mod_positions = [i for i, c in enumerate(pattern) if c != 'd']
    if not mod_positions:
        return -1
    
    length = len(pattern)
    if length == 0:
        return 0.0
    std = np.std(mod_positions)
    return std / length
###########################################################################

def compute_mod_block_count(pattern):
    """
    Returns the number of contiguous blocks of modified residues (non-'d')
    in the pattern.
    """
    pattern = str(pattern)
    blocks = re.findall(r'[^d]+', pattern)
    length = len(pattern)
    if length == 0:
        return 0.0

    return (len(blocks)) / length
############################################################################

def compute_mod_max_block_length(pattern):
    """
    Returns the length of the longest contiguous block of modified residues (non-'d')
    in the chemical pattern.
    """
    pattern = str(pattern)
    blocks = re.findall(r'[^d]+', pattern)
    return max((len(block) for block in blocks), default=0)
############################################################################


def compute_mod_char_entropy(pattern):
    """
    Computes the Shannon entropy of the non-'d' characters (modification types) in the pattern.
    Entropy is 0 if there are no modifications or if all are of the same type.
    """
    pattern = str(pattern)
    mod_chars = [c for c in pattern if c != 'd']
    if not mod_chars:
        return 0.0

    counts = Counter(mod_chars)
    total = sum(counts.values())
    probs = [count / total for count in counts.values()]
    entropy = -sum(p * np.log2(p) for p in probs)

    max_entropy = np.log2(len(probs)) if len(probs) > 1 else 1  # to avoid /0
    return entropy / max_entropy
############################################################################

def compute_dominant_mod_fraction(pattern):
    """
    Returns the relative frequency of the most common modified residue
    (non-'d') in the chemical pattern. Returns 0.0 if there are no modifications.
    """
    pattern = str(pattern)
    mod_chars = [c for c in pattern if c != 'd']
    if not mod_chars:
        return 0.0

    counts = Counter(mod_chars)
    dominant_count = counts.most_common(1)[0][1]
    return dominant_count / len(mod_chars)
############################################################################


def compute_mod_evenness(pattern):
    """
    Computes the normalized entropy of the distances between modified residues
    (non-'d') as a proxy for spatial evenness. Returns 0.0 if < 2 modifications.
    """
    pattern = str(pattern)
    positions = [i for i, c in enumerate(pattern) if c != 'd']
    
    if len(positions) < 2:
        return 0.0  # no spacing to measure

    spacings = np.diff(positions)
    counts = np.bincount(spacings)
    probs = counts[counts > 0] / counts.sum()
    entropy = -np.sum(probs * np.log2(probs))

    max_entropy = np.log2(len(spacings)) if len(spacings) > 1 else 1
    return entropy / max_entropy
############################################################################
def compute_mod_symmetry_score(pattern):
    """
    Computes a symmetry score of modification distribution around the center.
    Returns a value between 0 (fully asymmetric) and 1 (perfect symmetry).
    """
    pattern = str(pattern)
    total_length = len(pattern)
    if total_length == 0:
        return 1.0  # define empty pattern as trivially symmetric

    mid = total_length // 2
    left = pattern[:mid]
    right = pattern[-mid:]

    left_count = sum(1 for c in left if c != 'd')
    right_count = sum(1 for c in right if c != 'd')
    total = left_count + right_count

    if total == 0:
        return 1.0  # no modifications = trivially symmetric

    return 1 - abs(left_count - right_count) / total
############################################################################
def compute_mod_skew_index(pattern):
    """
    Computes skew index between the 5' and 3' thirds of the pattern.
    Returns a value between -1 (all in 3') and +1 (all in 5').
    """
    pattern = str(pattern)
    n = len(pattern)
    if n == 0:
        return 0.0

    third = n // 3
    five_prime = pattern[:third]
    three_prime = pattern[-third:]

    mod_5p = sum(1 for c in five_prime if c != 'd')
    mod_3p = sum(1 for c in three_prime if c != 'd')

    total = mod_5p + mod_3p
    if total == 0:
        return 0.0

    return (mod_5p - mod_3p) / total
############################################################################
def compute_mod_mean_gap(pattern):
    """
    Returns the mean gap between consecutive modified residues (non-'d').
    Returns -1 if fewer than 2 modifications are found.
    """
    pattern = str(pattern)
    positions = [i for i, c in enumerate(pattern) if c != 'd']
    
    if len(positions) < 2:
        return -1

    gaps = np.diff(positions)

    length = len(pattern)
    if length == 0:
        return 0        
    return (gaps.mean())/length
############################################################################
def compute_mod_local_density_max(pattern, window=5):
    """
    Returns the maximum number of modified residues (non-'d') found in any sliding window
    of length `window` across the given chemical pattern.
    """
    pattern = str(pattern)
    n = len(pattern)
    
    if n < window:
        return sum(1 for c in pattern if c != 'd')  # entire pattern

    max_count = 0
    for i in range(n - window + 1):
        window_seq = pattern[i:i+window]
        count = sum(1 for c in window_seq if c != 'd')
        if count > max_count:
            max_count = count

    return max_count / window
############################################################################

def compute_mod_in_core(pattern, core_fraction=0.4):
    """
    Returns the fraction of modified residues (non-'d') that lie within
    the central portion of the pattern (e.g., 30%-70%).
    Returns -1 if no modifications exist.
    """
    pattern = str(pattern)
    n = len(pattern)
    if n == 0:
        return -1

    core_start = int(n * (0.5 - core_fraction / 2))
    core_end = int(n * (0.5 + core_fraction / 2))

    total_mods = sum(1 for c in pattern if c != 'd')
    if total_mods == 0:
        return -1

    core_mods = sum(1 for c in pattern[core_start:core_end] if c != 'd')
    return core_mods / total_mods
#############################################################################
def compute_mod_longest_repeat_run(pattern):
    """
    Returns the length of the longest run of identical modified characters (non-'d')
    appearing consecutively in the pattern.
    """
    pattern = str(pattern)
    max_run = 0
    current_run = 0
    prev_char = None

    for c in pattern:
        if c == 'd':
            current_run = 0
            prev_char = None
            continue

        if c == prev_char:
            current_run += 1
        else:
            current_run = 1
            prev_char = c

        if current_run > max_run:
            max_run = current_run
    length = len(pattern)
    if length == 0:
        return 0  # define empty pattern as having no runs
    
    return max_run / length  # normalize by pattern length for consistency
############################################################################
def compute_mod_adjacent_pair_count(pattern):
    """
    Returns the number of adjacent identical modification pairs (e.g., 'LL', 'OO'),
    where both characters are the same and not 'd'.
    Each overlapping identical pair is counted once.
    """
    pattern = str(pattern)
    count = 0
    for i in range(len(pattern) - 1):
        if pattern[i] != 'd' and pattern[i] == pattern[i + 1]:
            count += 1
    
    length = len(pattern)
    if length == 0:
        return 0        
    return count/ length

############################################################################
def compute_mod_strong_repeat_group_count(pattern, min_run_length=3):
    """
    Returns the number of groups of identical, consecutive modification characters (non-'d'),
    where each group has length ≥ min_run_length (default 3).
    """
    pattern = str(pattern)
    count = 0
    i = 0
    n = len(pattern)

    while i < n - 1:
        if pattern[i] == 'd':
            i += 1
            continue

        run_length = 1
        while i + run_length < n and pattern[i] == pattern[i + run_length]:
            run_length += 1

        if run_length >= min_run_length:
            count += 1
            i += run_length
        else:
            i += 1

    length = len(pattern)
    if length == 0:
        return 0        
    return count/ length
############################################################################
