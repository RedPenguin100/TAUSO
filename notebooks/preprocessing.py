# Imports
import pandas as pd
import numpy as np
import tauso
from notebooks.consts import *
from notebooks.notebook_utils import log_correction, read_cached_gene_to_data
from tauso.util import get_antisense


def process_row(row,gene_to_data):
    gene_name = row[CANONICAL_GENE]
    if gene_name not in gene_to_data:
        return -1, 0, "", ""

    locus_info = gene_to_data[gene_name]
    pre_mrna = locus_info.full_mrna
    antisense = getattr(row, SEQUENCE)
    sense = get_antisense(antisense)
    idx = pre_mrna.find(sense)  # the index in the pre_mrna the sense is found

    return idx, len(antisense), sense, pre_mrna


def preprocess_aso_data(csv_path):
    """
    Loads raw ASO data, filters for human genes, removes controls,
    calculates log inhibition, and maps sequences to pre-mRNA.

    Returns:
        pd.DataFrame: Cleaned and enriched DataFrame.
    """

    # 1. Load Data
    if not csv_path.exists():
        raise FileNotFoundError(f"Data file not found at: {csv_path}")

    all_data = pd.read_csv(str(csv_path), low_memory=False)

    # 2. Basic Filtering & Log Correction
    df = all_data.dropna(subset=[INHIBITION]).copy()
    log_correction(df)

    # 3. Filter Human Only
    df = df[df[CELL_LINE_ORGANISM] == 'human'].copy()

    # 4. Filter Specific Genes
    genes = df[CANONICAL_GENE].unique().tolist()
    genes_to_remove = ['HBV', 'negative_control']
    genes_clean = [g for g in genes if g not in genes_to_remove]
    gene_to_data = read_cached_gene_to_data(genes_clean)
    df = df[df[CANONICAL_GENE].isin(genes_clean)].copy()

    # 5. Sequence Mapping (Mapping ASO to pre-mRNA)

    # Initialize new columns
    df[SENSE_SEQUENCE] = ""
    df[PRE_MRNA_SEQUENCE] = ""
    df[SENSE_START] = -1  # initialize to -1 to find errors
    df[SENSE_LENGTH] = 0

    results = df.apply(lambda row: process_row(row,gene_to_data), axis=1)

    df[[SENSE_START, SENSE_LENGTH, SENSE_SEQUENCE, PRE_MRNA_SEQUENCE]] = pd.DataFrame(results.tolist(), index=df.index)

    # 6. Final Cleanup - removing (idx == -1)
    valid_data = df[df[SENSE_START] != -1].copy()

    print(f"Preprocessing complete. Final valid rows: {len(valid_data)}")
    return valid_data