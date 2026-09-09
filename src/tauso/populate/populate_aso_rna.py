"""Geometry of the ASO:RNA duplex, averaged over each region of the gapmer."""

import logging

from ..data.consts import ASO_SEQUENCE, CHEMICAL_PATTERN
from ..features.aso_rna.aso_rna import FEATURE_NAMES, calculate_aso_rna
from ..pandas_utils import add_columns

logger = logging.getLogger(__name__)

PREFIX = "asorna"


def aso_rna_feature_name(quantity):
    """Stable column name for one observable and region."""
    if quantity not in FEATURE_NAMES:
        raise ValueError(f"unknown quantity {quantity!r}")
    return f"{PREFIX}_{quantity}"


def aso_rna_feature_names():
    return [aso_rna_feature_name(q) for q in FEATURE_NAMES]


def populate_aso_rna_features(df, cpus=1):
    """Add one column per observable, region and statistic."""
    missing = [c for c in (ASO_SEQUENCE, CHEMICAL_PATTERN) if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in DataFrame: {missing}")

    scored = calculate_aso_rna(df[ASO_SEQUENCE].to_numpy(), df[CHEMICAL_PATTERN].to_numpy())
    columns = {aso_rna_feature_name(q): scored[q] for q in FEATURE_NAMES}
    return add_columns(df, columns), list(columns)
