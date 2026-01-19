import numpy as np
import pandas as pd
from notebooks.consts import *
from notebooks.notebook_utils import log_correction, read_cached_gene_to_data
from tauso.new_model.utils import SMILES
from tauso.util import get_antisense


def get_unique_genes(df):
    """
    Extracts unique genes from a dataframe, excluding blacklist items.
    """
    genes = df[CANONICAL_GENE].unique()
    genes_to_remove = {'HBV', 'negative_control', np.nan}

    # Filter using list comprehension
    genes_clean = [g for g in genes if g not in genes_to_remove]
    return genes_clean


def get_filtered_human_data(df):
    """
    Filters the raw DataFrame for Human entries and removes rows with missing Inhibition data.
    Returns a copy to avoid SettingWithCopy warnings.
    """
    # Filter: Organism is Human AND Inhibition is not NaN
    condition = (df[CELL_LINE_ORGANISM] == 'human') & (df[INHIBITION].notna())
    return df[condition].copy()


def process_row(row, gene_to_data):
    gene_name = row[CANONICAL_GENE]
    if gene_name not in gene_to_data:
        return -1, 0, ""

    locus_info = gene_to_data[gene_name]
    pre_mrna = locus_info.full_mrna
    antisense = getattr(row, SEQUENCE)
    sense = get_antisense(antisense)
    idx = pre_mrna.find(sense)  # the index in the pre_mrna the sense is found

    return idx, len(antisense), sense


def preprocess_aso_data(csv_path, include_smiles: bool = False):
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

    # 3. Filter Specific Genes (using the clean DF)
    genes_clean = get_unique_genes(df)
    gene_to_data = read_cached_gene_to_data(genes_clean)
    df = df[df[CANONICAL_GENE].isin(genes_clean)].copy()

    # 4. Optional: Drop SMILES
    if not include_smiles and SMILES in df.columns:
        df.drop(columns=[SMILES], inplace=True)

    # 5. Sequence Mapping
    # Initialize columns
    df[SENSE_SEQUENCE] = ""
    df[SENSE_START] = -1
    df[SENSE_LENGTH] = 0

    results = df.apply(lambda row: process_row(row, gene_to_data), axis=1)
    df[[SENSE_START, SENSE_LENGTH, SENSE_SEQUENCE]] = pd.DataFrame(results.tolist(), index=df.index)

    # 6. Final Cleanup
    valid_data = df[df[SENSE_START] != -1].copy()

    print(f"Preprocessing complete. Final valid rows: {len(valid_data)}")
    return valid_data