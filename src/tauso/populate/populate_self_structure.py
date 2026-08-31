import logging

from ..data.consts import ASO_SEQUENCE, CHEMICAL_PATTERN
from ..features.self_structure.self_structure import FEATURE_NAMES, calculate_self_structure

logger = logging.getLogger(__name__)

# The weight set the nearest-neighbour tables were fitted with. It is carried in the column
# name because a refit changes the numbers under a feature that otherwise looks the same.
WEIGHT_SET = "md3"


def self_structure_feature_name(quantity, weight_set=WEIGHT_SET):
    """Stable column name for one self-structure quantity."""
    if quantity not in FEATURE_NAMES:
        raise ValueError(f"unknown quantity {quantity!r}; expected one of {', '.join(FEATURE_NAMES)}")
    return f"{weight_set}_{quantity}"


def self_structure_feature_names(weight_set=WEIGHT_SET):
    return [self_structure_feature_name(q, weight_set) for q in FEATURE_NAMES]


def populate_self_structure_features(df, weight_set=WEIGHT_SET):
    """Add one column per self-structure quantity.

    Scored from the oligo alone -- its sequence and its chemistry -- so no target, locus or
    fold is involved.
    """
    missing = [c for c in (ASO_SEQUENCE, CHEMICAL_PATTERN) if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in DataFrame: {missing}")

    scored = calculate_self_structure(df[ASO_SEQUENCE].to_numpy(), df[CHEMICAL_PATTERN].to_numpy())
    names = []
    for quantity in FEATURE_NAMES:
        name = self_structure_feature_name(quantity, weight_set)
        df[name] = scored[quantity]
        names.append(name)
    return df, names
