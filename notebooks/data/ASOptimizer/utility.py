import re

import numpy as np
import pandas as pd

from notebooks.data.OligoAI.parse_chemistry import assign_backbone
from notebooks.data.utility import standardize_cell_line
from tauso.data.consts import (
    CANONICAL_GENE,
    SEQUENCE,
    LINKAGE_LOCATION,
    TRANSFECTION,
)


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


def transfection_to_oligo(df):
    # Ensure the Transfection column is string type and handle NaNs
    col = df[TRANSFECTION].astype(str).str.lower()

    # Define conditions
    # 1. Gymnosis: matches "free uptake" or "uptake"
    is_gymnosis = col.str.contains("uptake", na=False)

    # 2. Electroporation: exact or partial match
    is_electroporation = col.str.contains("electroporation", na=False)

    # 3. Lipofection: contains "lipo" OR "cytofectin"
    is_lipofection = col.str.contains("lipo|Cytofectin", na=False)

    # Apply logic using np.select (cleanest way for multiple conditions)
    condlist = [is_gymnosis, is_electroporation, is_lipofection]
    choicelist = ["Gymnosis", "Electroporation", "Lipofection"]

    df["transfection_method"] = np.select(condlist, choicelist, default="Other")

    return df


def standardize_asoptimizer_data(df: pd.DataFrame) -> pd.DataFrame:
    result = df.copy()

    # HBV is not handled well in this version of the code. # TODO: fix HBV
    result = result[result[CANONICAL_GENE] != "HBV"].reset_index(drop=True)

    # Convert ASOptimizer transfection to OligoAI transfection
    result = transfection_to_oligo(result)

    result["backbone_mods"] = result.apply(
        lambda row: transform_linkage_to_oligo(row[LINKAGE_LOCATION], len(row[SEQUENCE])), axis=1
    )

    result = assign_backbone(result)

    result = standardize_cell_line(result)
    return result