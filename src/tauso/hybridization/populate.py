import pandas as pd
import multiprocessing
from pandarallel import pandarallel

from tauso.hybridization.hybridization_features import calc_methylcytosines, calculate_lna, calculate_cet, \
    get_exp_psrna_hybridization, get_exp_psrna_hybridization_diff, calculate_dna, get_exp_dna_rna_hybridization
from tauso.hybridization.md_weights import get_2moe_md_diff, get_psdna_rna_md_total
from tauso.new_model.consts_dataframe import *

# Logic can be a lambda (row-wise) or None (for vectorized post-processing)
HYBR_FEATURE_TO_CALCULATION = {
    'MOE_DIFF_37_MD_GB_HYBR': lambda row: get_2moe_md_diff(
        row[SEQUENCE], row[CHEMICAL_PATTERN], row[MODIFICATION], simul_type='gb'
    ),
    'MOE_DIFF_37_MD_PB_HYBR': lambda row: get_2moe_md_diff(
        row[SEQUENCE], row[CHEMICAL_PATTERN], row[MODIFICATION], simul_type='pb'
    ),
    'PSDNA_RNA_MD_37_GB_TOTAL_HYBR': lambda row: get_psdna_rna_md_total(
        row[SEQUENCE], row[MODIFICATION], simul_type='gb'
    ),
    'PSDNA_RNA_MD_37_PB_TOTAL_HYBR': lambda row: get_psdna_rna_md_total(
        row[SEQUENCE], row[MODIFICATION], simul_type='pb'
    ),
    'METHYL_CYTOSINES': lambda row: calc_methylcytosines(
        row[SEQUENCE], row[CHEMICAL_PATTERN]
    ),
    'LNA_DIFF_37_HYBR': lambda row: calculate_lna(
        row[SEQUENCE], row[CHEMICAL_PATTERN]
    ),
    'CET_DIFF_37_HYBR': lambda row: calculate_cet(
        row[SEQUENCE], row[CHEMICAL_PATTERN]
    ),
    'TOTAL_PSDNA_HYBR': lambda row: get_exp_psrna_hybridization(row[SEQUENCE]) / 1000,
    'PSDNA_DIFF_37_HYBR': lambda row: get_exp_psrna_hybridization_diff(row[SEQUENCE]) / 1000,
    'TOTAL_DNA_HYBR': lambda row: calculate_dna(row[SEQUENCE]),
    'TOTAL_DNA_RNA_HYBR': lambda row: get_exp_dna_rna_hybridization(row[SEQUENCE]),

    # Derived features (Calculated via vectorization later)
    'DNA_HYBR_DIFF': None
}


def populate_hybridization(df, n_cores=1, features_to_run=None):
    """
    Populates hybridization features using the FEATURE_TO_CALCULATION registry.
    """
    all_data = df.copy()

    # Default to all features in registry if none specified
    if features_to_run is None:
        features_to_run = list(HYBR_FEATURE_TO_CALCULATION.keys())

    # 1. Setup Parallel Engine
    nb_workers = n_cores if n_cores is not None else multiprocessing.cpu_count()
    if nb_workers > 1:
        pandarallel.initialize(nb_workers=nb_workers, progress_bar=False)
        apply_func = all_data.parallel_apply
    else:
        apply_func = all_data.apply

    # 2. Run Row-wise Calculations (Parallelized)
    for feature in features_to_run:
        logic = HYBR_FEATURE_TO_CALCULATION.get(feature)

        # Only run if logic is a function/lambda
        if callable(logic):
            print(f"Calculating feature: {feature}...")
            all_data[feature] = apply_func(logic, axis=1)

    # 3. Run Vectorized Calculations (Fast Post-processing)
    # This avoids redundant function calls for simple arithmetic
    if 'DNA_HYBR_DIFF' in features_to_run:
        print("Calculating vectorized feature: DNA_HYBR_DIFF...")
        all_data['DNA_HYBR_DIFF'] = all_data['TOTAL_DNA_HYBR'] - all_data['TOTAL_DNA_RNA_HYBR']

    return all_data