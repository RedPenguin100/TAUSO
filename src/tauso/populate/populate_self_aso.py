"""Features of the oligo on its own: what it folds into, and what it does with a second copy.

Everything here is scored from the oligo's sequence and chemistry alone. No target, no locus,
no fold of anything but the oligo, which is what separates this family from the sequence
features it used to sit among.

Three searches contribute (see `features.self_aso.self_aso`), plus the oligo folded with RNA
parameters, which is kept under its original column name so existing matrices still line up.
"""

import logging

from ..data.consts import ASO_SEQUENCE, CHEMICAL_PATTERN
from ..features.self_aso.self_aso import FEATURE_NAMES, calculate_self_aso
from ..features.sequence.seq_features import internal_fold_rna
from ..parallel_utils import make_apply_fn

logger = logging.getLogger(__name__)

PREFIX = "selfaso"

RNA_FOLD_FEATURE = "seq_internal_fold_rna"
"""The oligo folded with RNA (Turner) parameters.

Physically the wrong alphabet for a gapmer, and deliberately so: it finds structure in 38% of
oligos where the chemistry-aware search finds a hairpin in 6%, and that sensitivity is what
makes it useful. Kept under its original name because renaming it would break every stored
matrix and the interaction built on it.
"""


def self_aso_feature_name(quantity):
    """Stable column name for one self-structure quantity."""
    if quantity not in FEATURE_NAMES:
        raise ValueError(f"unknown quantity {quantity!r}; expected one of {', '.join(FEATURE_NAMES)}")
    return f"{PREFIX}_{quantity}"


def self_aso_feature_names():
    return [self_aso_feature_name(q) for q in FEATURE_NAMES] + [RNA_FOLD_FEATURE]


def populate_self_aso_features(df, cpus=1):
    """Add one column per self-structure quantity, plus the RNA fold of the oligo."""
    missing = [c for c in (ASO_SEQUENCE, CHEMICAL_PATTERN) if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in DataFrame: {missing}")

    scored = calculate_self_aso(df[ASO_SEQUENCE].to_numpy(), df[CHEMICAL_PATTERN].to_numpy())
    names = []
    for quantity in FEATURE_NAMES:
        name = self_aso_feature_name(quantity)
        df[name] = scored[quantity]
        names.append(name)

    df[RNA_FOLD_FEATURE] = make_apply_fn(df[ASO_SEQUENCE], n_jobs=cpus)(internal_fold_rna)
    names.append(RNA_FOLD_FEATURE)
    return df, names
