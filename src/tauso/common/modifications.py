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