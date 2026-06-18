import re


def compute_mod_sugar_block_count(pattern):
    """
    Returns the number of contiguous blocks of modified residues (non-'d')
    in the pattern.
    """
    pattern = str(pattern)
    return len(re.findall(r"[^d]+", pattern))


def compute_mod_sugar_max_block_length(pattern):
    """
    Returns the length of the longest contiguous block of modified residues (non-'d')
    in the chemical pattern.
    """
    pattern = str(pattern)
    blocks = re.findall(r"[^d]+", pattern)
    return max((len(block) for block in blocks), default=0)
