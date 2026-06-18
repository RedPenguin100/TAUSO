import logging
import multiprocessing

logger = logging.getLogger(__name__)

from ..data.consts import CHEMICAL_PATTERN
from ..features.sequence_modification.mod_features import (
    compute_mod_sugar_block_count,
    compute_mod_sugar_max_block_length,
)
from ..parallel_utils import make_apply_fn

MODIFICATION_FEATURE_TO_CALCULATION = {
    "mod_sugar_block_count": lambda row: compute_mod_sugar_block_count(row[CHEMICAL_PATTERN]),
    "mod_sugar_max_block_length": lambda row: compute_mod_sugar_max_block_length(row[CHEMICAL_PATTERN]),
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
