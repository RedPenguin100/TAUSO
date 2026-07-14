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


def transform_linkage_to_oligo(pattern_string, seq_length):
    """
    Transforms linkage shorthand into a string representation of a list.
    Length of PS/PO list is seq_length - 1, followed by one <PAD>.
    """
    # 1. Determine number of linkage slots (n-1)
    num_linkages = seq_length - 1

    # 2. Handle the 'else' case: everything is PS
    if str(pattern_string).strip().lower() == "else":
        result = ["PS"] * num_linkages
    else:
        # 3. Extract PS indices from string (e.g., "0?4?15" -> {0, 4, 15})
        ps_indices = set(map(int, re.findall(r"\d+", str(pattern_string))))

        # 4. Build the list of PS/PO
        result = []
        for i in range(num_linkages):
            if i in ps_indices:
                result.append("PS")
            else:
                result.append("PO")

    # 5. Always add the PAD at the end
    result.append("<PAD>")

    # 6. Return as the exact string format requested
    return str(result)


def transform_pattern_to_oligo(sequence_string, lna_as_cet=True):
    """
    Transforms shorthand to a STRING representation of a list.
    Input: "MMM" -> Output: "['MOE', 'MOE', 'MOE']"
    """
    mapping = {"M": "MOE", "d": "DNA", "C": "CET"}
    if lna_as_cet:
        mapping["L"] = "CET"

    # Create the list first
    result_list = [mapping.get(char, char) for char in sequence_string]

    # Convert the list object into a literal string
    return str(result_list)
