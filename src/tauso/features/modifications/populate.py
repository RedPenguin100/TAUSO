from tauso.features.mod_features import (
    compute_mod_fraction, compute_mod_type_count, compute_mod_5prime_run,
    compute_mod_3prime_run, compute_mod_min_distance_to_5prime,
    compute_mod_min_distance_to_3prime, compute_mod_pos_std,
    compute_mod_block_count, compute_mod_max_block_length,
    compute_mod_char_entropy, compute_dominant_mod_fraction,
    compute_mod_evenness, compute_mod_symmetry_score,
    compute_mod_skew_index, compute_mod_mean_gap,
    compute_mod_local_density_max, compute_mod_in_core,
    compute_mod_longest_repeat_run, compute_mod_adjacent_pair_count,
    compute_mod_strong_repeat_group_count
)
from tauso.new_model.consts_dataframe import CHEMICAL_PATTERN

import pandas as pd
import multiprocessing
from pandarallel import pandarallel

MODIFICATION_FEATURE_TO_CALCULATION = {
    'Modification_fraction': lambda row: compute_mod_fraction(row[CHEMICAL_PATTERN]),
    'Modification_type_count': lambda row: compute_mod_type_count(row[CHEMICAL_PATTERN]),
    'Modification_5prime_run': lambda row: compute_mod_5prime_run(row[CHEMICAL_PATTERN]),
    'Modification_3prime_run': lambda row: compute_mod_3prime_run(row[CHEMICAL_PATTERN]),
    'Modification_min_distance_to_5prime': lambda row: compute_mod_min_distance_to_5prime(row[CHEMICAL_PATTERN]),
    'Modification_min_distance_to_3prime': lambda row: compute_mod_min_distance_to_3prime(row[CHEMICAL_PATTERN]),
    'Modification_pos_std': lambda row: compute_mod_pos_std(row[CHEMICAL_PATTERN]),
    'Modification_block_count': lambda row: compute_mod_block_count(row[CHEMICAL_PATTERN]),
    'Modification_max_block_length': lambda row: compute_mod_max_block_length(row[CHEMICAL_PATTERN]),
    'Modification_char_entropy': lambda row: compute_mod_char_entropy(row[CHEMICAL_PATTERN]),
    'Modification_dominant_mod_fraction': lambda row: compute_dominant_mod_fraction(row[CHEMICAL_PATTERN]),
    'Modification_evenness': lambda row: compute_mod_evenness(row[CHEMICAL_PATTERN]),
    'Modification_symmetry_score': lambda row: compute_mod_symmetry_score(row[CHEMICAL_PATTERN]),
    'Modification_skew_index': lambda row: compute_mod_skew_index(row[CHEMICAL_PATTERN]),
    'Modification_mean_gap': lambda row: compute_mod_mean_gap(row[CHEMICAL_PATTERN]),
    'Modification_local_density_max': lambda row: compute_mod_local_density_max(row[CHEMICAL_PATTERN]),
    'Modification_in_core': lambda row: compute_mod_in_core(row[CHEMICAL_PATTERN]),
    'Modification_longest_repeat_run': lambda row: compute_mod_longest_repeat_run(row[CHEMICAL_PATTERN]),
    'Modification_adjacent_pair_count': lambda row: compute_mod_adjacent_pair_count(row[CHEMICAL_PATTERN]),
    'Modification_strong_repeat_group_count': lambda row: compute_mod_strong_repeat_group_count(row[CHEMICAL_PATTERN]),
}


def populate_modifications(df, n_cores=None):
    """
    Populates modification-based features for ASO chemical patterns.

    Args:
        df (pd.DataFrame): Input dataframe.
        n_cores (int/None): Number of cores to use. None uses all available.
    """
    all_data = df.copy()

    if n_cores is None:
        nb_workers = multiprocessing.cpu_count()
    else:
        nb_workers = n_cores

    # Initialize parallel logic or fallback to standard pandas
    if nb_workers > 1:
        # initialize is called once; it won't re-init if already active in most environments
        pandarallel.initialize(nb_workers=nb_workers)
        apply_func = all_data.parallel_apply
    else:
        apply_func = all_data.apply

    # 3. Execution Loop
    for feature in MODIFICATION_FEATURE_TO_CALCULATION.keys():
        logic = MODIFICATION_FEATURE_TO_CALCULATION.get(feature)

        if callable(logic):
            print(f"Calculating: {feature}")
            # Note: ensure CHEMICAL_PATTERN constant is accessible here
            all_data[feature] = apply_func(logic, axis=1)
        else:
            print(f"Warning: {feature} logic not found or not callable. Skipping.")

    return all_data