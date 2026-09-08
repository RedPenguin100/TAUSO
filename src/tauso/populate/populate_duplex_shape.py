"""Duplex geometry along the oligo, by region.

Scored from the oligo's sequence and chemistry alone -- no target and no folding. See
`features.duplex_shape.duplex_shape` for the lookup and how the regions are drawn.
"""

import logging

from ..data.consts import ASO_SEQUENCE, CHEMICAL_PATTERN
from ..features.duplex_shape.duplex_shape import QUANTITIES, calculate_duplex_shape
from ..pandas_utils import add_columns

logger = logging.getLogger(__name__)

PREFIX = "shape"


def duplex_shape_feature_name(quantity):
    """Stable column name for one duplex-shape quantity."""
    if quantity not in QUANTITIES:
        raise ValueError(f"unknown quantity {quantity!r}; expected one of {', '.join(QUANTITIES)}")
    return f"{PREFIX}_{quantity}"


def duplex_shape_feature_names():
    return [duplex_shape_feature_name(q) for q in QUANTITIES]


def populate_duplex_shape_features(df):
    """Add one column per duplex-shape quantity."""
    missing = [c for c in (ASO_SEQUENCE, CHEMICAL_PATTERN) if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in DataFrame: {missing}")

    scored = calculate_duplex_shape(df[ASO_SEQUENCE].to_numpy(), df[CHEMICAL_PATTERN].to_numpy())
    columns = {duplex_shape_feature_name(q): scored[q] for q in QUANTITIES}
    return add_columns(df, columns), list(columns)
