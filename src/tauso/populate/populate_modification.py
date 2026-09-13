import logging

logger = logging.getLogger(__name__)

from ..data.consts import CHEMICAL_PATTERN
from ..features.sequence_modification.mod_features import (
    compute_mod_sugar_block_count,
    compute_mod_sugar_max_block_length,
)
from ..pandas_utils import add_columns

MODIFICATION_FEATURE_TO_CALCULATION = {
    "mod_sugar_block_count": compute_mod_sugar_block_count,
    "mod_sugar_max_block_length": compute_mod_sugar_max_block_length,
}
"""Each reads the chemical pattern alone, so the features are mapped over that column."""


def populate_modifications(df, features_to_run=None):
    """
    Populates modification-based features for ASO chemical patterns.

    Args:
        df (pd.DataFrame): Input dataframe.
        features_to_run (list/None): Features to calculate. None runs the whole registry.
    """
    if features_to_run is None:
        features_to_run = list(MODIFICATION_FEATURE_TO_CALCULATION.keys())

    # 3. Execution Loop
    computed = {}
    for feature in features_to_run:
        logic = MODIFICATION_FEATURE_TO_CALCULATION.get(feature)

        if callable(logic):
            logger.debug("Calculating: %s", feature)
            computed[feature] = df[CHEMICAL_PATTERN].map(logic)
        else:
            logger.warning("Feature '%s' logic not found or not callable. Skipping.", feature)

    return add_columns(df, computed), features_to_run
