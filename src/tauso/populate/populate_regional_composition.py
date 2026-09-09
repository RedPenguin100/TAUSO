"""Regional dinucleotide composition of the oligo, relative to the oligo as a whole."""

import logging

from ..data.consts import ASO_SEQUENCE, CHEMICAL_PATTERN
from ..features.regional_composition.regional_composition import (
    FEATURE_NAMES,
    calculate_regional_composition,
)
from ..pandas_utils import add_columns

logger = logging.getLogger(__name__)

PREFIX = "dinuc"


def regional_composition_feature_name(quantity):
    """Stable column name for one region/dinucleotide pair."""
    if quantity not in FEATURE_NAMES:
        raise ValueError(f"unknown quantity {quantity!r}")
    return f"{PREFIX}_{quantity}"


def regional_composition_feature_names():
    return [regional_composition_feature_name(q) for q in FEATURE_NAMES]


def populate_regional_composition_features(df, cpus=1):
    """Add one column per region and dinucleotide."""
    missing = [c for c in (ASO_SEQUENCE, CHEMICAL_PATTERN) if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in DataFrame: {missing}")

    scored = calculate_regional_composition(df[ASO_SEQUENCE].to_numpy(), df[CHEMICAL_PATTERN].to_numpy())
    columns = {regional_composition_feature_name(q): scored[q] for q in FEATURE_NAMES}
    return add_columns(df, columns), list(columns)
