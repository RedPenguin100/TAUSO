import re

import numpy as np
import pandas as pd
from notebooks.data.ASOptimizer.consts import LINKAGE_LOCATION_ASOPT, SMILES_ASOPT, TRANSFECTION_ASOPT
from notebooks.data.OligoAI.parse_chemistry import assign_backbone
from notebooks.data.utility import standardize_cell_line
from notebooks.notebook_utils import get_unique_genes, log_correction
from notebooks.preprocessing import process_row

from tauso.data.consts import (
    ASO_SEQUENCE,
    CANONICAL_GENE_NAME,
    CELL_LINE,
    CELL_LINE_DEPMAP,
    CELL_LINE_DEPMAP_PROXY,
    CELL_LINE_ORGANISM,
    CHEMICAL_PATTERN,
    INHIBITION_PERCENT,
    STRUCTURE_SENSE_LENGTH,
    STRUCTURE_SENSE_SEQUENCE,
    STRUCTURE_SENSE_START,
    resolve_depmap_id,
    resolve_depmap_proxy,
)
from tauso.genome.read_human_genome import get_locus_to_data_dict


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
    col = df[TRANSFECTION_ASOPT].astype(str).str.lower()

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
    result = result[result[CANONICAL_GENE_NAME] != "HBV"].reset_index(drop=True)

    # Convert ASOptimizer transfection to OligoAI transfection
    result = transfection_to_oligo(result)

    result["backbone_mods"] = result.apply(
        lambda row: transform_linkage_to_oligo(row[LINKAGE_LOCATION_ASOPT], len(row[ASO_SEQUENCE])), axis=1
    )

    result = assign_backbone(result)

    result = standardize_cell_line(result)
    return result


def get_filtered_human_data(df):
    """
    Filters the raw DataFrame for Human entries and removes rows with missing Inhibition data.
    Returns a copy to avoid SettingWithCopy warnings.
    """
    # Filter: Organism is Human AND Inhibition is not NaN
    condition = (df[CELL_LINE_ORGANISM] == "human") & (df[INHIBITION_PERCENT].notna())
    return df[condition].copy()


def preprocess_asoptimizer_data(csv_path, include_smiles: bool = False):
    """
    Loads raw ASO data, cleans it, calculates log inhibition, and maps sequences.
    """
    if not csv_path.exists():
        raise FileNotFoundError(f"Data file not found at: {csv_path}")

    # 1. Load Data
    all_data = pd.read_csv(str(csv_path), low_memory=False)

    # 2. Filter & Log Correct (Abstracted filtering logic)
    df = get_filtered_human_data(all_data)
    log_correction(df)

    df = df[df[CHEMICAL_PATTERN] != "(S)-cEt/5-methylcytosines/deoxy"].copy()
    df = df[df[CHEMICAL_PATTERN] != "LNA/deoxy"].copy()

    # 3. Filter Specific Genes (using the clean DF)
    genes_clean = get_unique_genes(df)
    gene_to_data = get_locus_to_data_dict(gene_subset=genes_clean)
    df = df[df[CANONICAL_GENE_NAME].isin(genes_clean)].copy()

    # 4. Optional: Drop SMILES
    if not include_smiles and SMILES_ASOPT in df.columns:
        df.drop(columns=[SMILES_ASOPT], inplace=True)

    # 5. Sequence Mapping
    # Initialize columns
    df[STRUCTURE_SENSE_SEQUENCE] = ""
    df[STRUCTURE_SENSE_START] = -1
    df[STRUCTURE_SENSE_LENGTH] = 0

    results = df.apply(lambda row: process_row(row, gene_to_data), axis=1)
    df[[STRUCTURE_SENSE_START, STRUCTURE_SENSE_LENGTH, STRUCTURE_SENSE_SEQUENCE]] = pd.DataFrame(results.tolist(), index=df.index)

    # 6. Final Cleanup
    valid_data = df[df[STRUCTURE_SENSE_START] != -1].copy()

    # 7. Add standard cell line column
    if not CELL_LINE_DEPMAP_PROXY in valid_data.columns:
        print("Adding cell line depmap name proxy")
        valid_data[CELL_LINE_DEPMAP_PROXY] = valid_data[CELL_LINE].map(resolve_depmap_proxy)

    # 8. Add Depmap cell lines
    if not CELL_LINE_DEPMAP in valid_data.columns:
        print("Adding a depmap column")
        valid_data[CELL_LINE_DEPMAP] = valid_data[CELL_LINE].map(resolve_depmap_id)

    print(f"Preprocessing complete. Final valid rows: {len(valid_data)}")
    return valid_data