import ast
import logging
import re
from functools import lru_cache

import pandas as pd

from tauso.data.consts import (
    BACKBONE_MODS,
    CHEMICAL_PATTERN,
    MIXMER_MODIFICATION,
    MODIFICATION_STRING,
    PS_PATTERN,
    SUGAR_MODS,
)

logger = logging.getLogger(__name__)

# One character per residue. Every sugar the source data can contain must appear here: the
# pattern is positional, so a residue with no character silently shifts every residue after it
# and leaves the pattern shorter than the oligo.
SUGAR_CHAR = {
    "MOE": "M",
    "CET": "C",
    "LNA": "L",
    "DNA": "d",
    "OME": "o",  # 2'-O-methyl
    "F": "f",  # 2'-fluoro
}

# The high-affinity sugar a gapmer is named for. An oligo carrying any other sugar, or more
# than one of these, is a mixmer and is excluded upstream.
GAPMER_LABEL = {
    "MOE": "MOE/5-methylcytosines/deoxy",
    "CET": "cEt/5-methylcytosines/deoxy",
    "LNA": "LNA/5-methylcytosines/deoxy",
}


# Cache will help store the most common modifications
@lru_cache(maxsize=100)
def _process_chemistry(mod_str):
    if not isinstance(mod_str, str) or not mod_str.startswith("["):
        return None, None

    try:
        mod_list = ast.literal_eval(mod_str)
        mods = [m.upper() for m in mod_list]

        unknown = sorted({m for m in mods if m not in SUGAR_CHAR})
        if unknown:
            logger.warning("Unknown sugar %s in %r; leaving the chemistry unparsed.", unknown, mod_str)
            return None, None

        # 1. Generate the CHEMICAL_PATTERN (e.g., MMMMMddddddddddMMMMM)
        pattern = "".join(SUGAR_CHAR[m] for m in mods)

        # 2. Determine the MODIFICATION_STRING label
        sugars = set(mods) - {"DNA"}
        gapmer_sugars = sugars & set(GAPMER_LABEL)
        if len(sugars) > 1 or sugars - gapmer_sugars:
            label = MIXMER_MODIFICATION
        elif gapmer_sugars:
            label = GAPMER_LABEL[next(iter(gapmer_sugars))]
        else:
            label = "DNA"

        return pattern, label
    except (ValueError, SyntaxError):
        return None, None


# Cache will help store the most common modifications
@lru_cache(maxsize=100)
def _parse_backbone(backbone_str):
    if not isinstance(backbone_str, str) or not backbone_str.startswith("["):
        return None

    try:
        raw_list = ast.literal_eval(backbone_str)
        mapping = {"PS": "*", "PO": "d", "<PAD>": ""}
        return "".join(mapping.get(item, "") for item in raw_list)
    except (ValueError, SyntaxError):
        return None


def assign_sugar(df: pd.DataFrame) -> pd.DataFrame:
    """Parses sugar chemistry strings and populates chemical pattern and modification labels."""
    if SUGAR_MODS not in df.columns:
        raise ValueError(f"column {SUGAR_MODS} is missing from the data, please populate it")

    # Assuming _process_chemistry is defined locally or imported
    results = df[SUGAR_MODS].apply(_process_chemistry)
    df[[CHEMICAL_PATTERN, MODIFICATION_STRING]] = pd.DataFrame(results.tolist(), index=df.index)
    return df


def assign_backbone(df: pd.DataFrame) -> pd.DataFrame:
    """Parses backbone chemistry strings and populates PS pattern labels."""
    if BACKBONE_MODS not in df.columns:
        raise ValueError(f"column {BACKBONE_MODS} is missing from the data, please populate it")

    # Assuming _parse_backbone is defined locally or imported
    df[PS_PATTERN] = df[BACKBONE_MODS].apply(_parse_backbone)
    return df


def assign_chemistry(df: pd.DataFrame) -> pd.DataFrame:
    """
    Parses sugar and backbone chemistry strings in a dataframe and populates new columns
    with their respective patterns and modification labels.
    """
    out_df = df.copy()

    out_df = assign_sugar(out_df)
    out_df = assign_backbone(out_df)

    return out_df


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
