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
