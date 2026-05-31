import logging
import multiprocessing

logger = logging.getLogger(__name__)

from ..data.consts import CHEMICAL_PATTERN
from ..features.sequence_modification.mod_features import (
    compute_mod_sugar_dominant_mod_fraction,
    compute_mod_sugar_3prime_run,
    compute_mod_sugar_5prime_run,
    compute_mod_sugar_adjacent_pair_count,
    compute_mod_sugar_block_count,
    compute_mod_sugar_char_entropy,
    compute_mod_sugar_evenness,
    compute_mod_sugar_fraction,
    compute_mod_sugar_in_core,
    compute_mod_sugar_local_density_max,
    compute_mod_sugar_longest_repeat_run,
    compute_mod_sugar_max_block_length,
    compute_mod_sugar_mean_gap,
    compute_mod_sugar_min_distance_to_3prime,
    compute_mod_sugar_min_distance_to_5prime,
    compute_mod_sugar_pos_std,
    compute_mod_sugar_skew_index,
    compute_mod_sugar_strong_repeat_group_count,
    compute_mod_sugar_symmetry_score,
)
from ..parallel_utils import make_apply_fn

MODIFICATION_FEATURE_TO_CALCULATION = {
    "mod_sugar_fraction": lambda row: compute_mod_sugar_fraction(row[CHEMICAL_PATTERN]),
    "mod_sugar_5prime_run": lambda row: compute_mod_sugar_5prime_run(row[CHEMICAL_PATTERN]),
    "mod_sugar_3prime_run": lambda row: compute_mod_sugar_3prime_run(row[CHEMICAL_PATTERN]),
    "mod_sugar_min_distance_to_5prime": lambda row: compute_mod_sugar_min_distance_to_5prime(row[CHEMICAL_PATTERN]),
    "mod_sugar_min_distance_to_3prime": lambda row: compute_mod_sugar_min_distance_to_3prime(row[CHEMICAL_PATTERN]),
    "mod_sugar_pos_std": lambda row: compute_mod_sugar_pos_std(row[CHEMICAL_PATTERN]),
    "mod_sugar_block_count": lambda row: compute_mod_sugar_block_count(row[CHEMICAL_PATTERN]),
    "mod_sugar_max_block_length": lambda row: compute_mod_sugar_max_block_length(row[CHEMICAL_PATTERN]),
    "mod_sugar_char_entropy": lambda row: compute_mod_sugar_char_entropy(row[CHEMICAL_PATTERN]),
    "mod_sugar_dominant_mod_fraction": lambda row: compute_mod_sugar_dominant_mod_fraction(row[CHEMICAL_PATTERN]),
    "mod_sugar_evenness": lambda row: compute_mod_sugar_evenness(row[CHEMICAL_PATTERN]),
    "mod_sugar_symmetry_score": lambda row: compute_mod_sugar_symmetry_score(row[CHEMICAL_PATTERN]),
    "mod_sugar_skew_index": lambda row: compute_mod_sugar_skew_index(row[CHEMICAL_PATTERN]),
    "mod_sugar_mean_gap": lambda row: compute_mod_sugar_mean_gap(row[CHEMICAL_PATTERN]),
    "mod_sugar_local_density_max": lambda row: compute_mod_sugar_local_density_max(row[CHEMICAL_PATTERN]),
    "mod_sugar_in_core": lambda row: compute_mod_sugar_in_core(row[CHEMICAL_PATTERN]),
    "mod_sugar_longest_repeat_run": lambda row: compute_mod_sugar_longest_repeat_run(row[CHEMICAL_PATTERN]),
    "mod_sugar_adjacent_pair_count": lambda row: compute_mod_sugar_adjacent_pair_count(row[CHEMICAL_PATTERN]),
    "mod_sugar_strong_repeat_group_count": lambda row: compute_mod_sugar_strong_repeat_group_count(row[CHEMICAL_PATTERN]),
}


def populate_modifications(df, n_cores=None, features_to_run=None):  # Added features_to_run
    """
    Populates modification-based features for ASO chemical patterns.

    Args:
        df (pd.DataFrame): Input dataframe.
        n_cores (int/None): Number of cores to use. None uses all available.
    """
    all_data = df.copy()

    nb_workers = multiprocessing.cpu_count() if n_cores is None else n_cores
    apply_func = make_apply_fn(all_data, n_jobs=nb_workers)

    if features_to_run is None:
        features_to_run = list(MODIFICATION_FEATURE_TO_CALCULATION.keys())

    # 3. Execution Loop
    for feature in features_to_run:
        logic = MODIFICATION_FEATURE_TO_CALCULATION.get(feature)

        if callable(logic):
            logger.debug("Calculating: %s", feature)
            all_data[feature] = apply_func(logic, axis=1)
        else:
            logger.warning("Feature '%s' logic not found or not callable. Skipping.", feature)

    return all_data, features_to_run
