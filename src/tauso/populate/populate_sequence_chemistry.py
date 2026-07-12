from typing import Iterable, Optional, Tuple

import pandas as pd

from ..data.consts import ASO_SEQUENCE, CHEMICAL_PATTERN
from ..features.sequence.seq_chemistry import gap_gc_content, wing_gap_gc_delta
from .feature_runner import FeatureSpec, compute_features

FEATURE_SPECS: list[FeatureSpec] = [
    ("mod_sugar_gap_gc_content", gap_gc_content),
    ("mod_sugar_wing_gap_gc_delta", wing_gap_gc_delta),
]


def populate_sequence_chemistry_features(
    df,
    features: Optional[Iterable[str]] = None,
    cpus: int = 1,
) -> Tuple:
    # Bundle the target columns into a single Series so features apply element-wise,
    # avoiding a per-row Series copy.
    targets = pd.Series(list(zip(df[ASO_SEQUENCE], df[CHEMICAL_PATTERN])), index=df.index, dtype=object)
    return compute_features(
        df,
        FEATURE_SPECS,
        targets,
        lambda apply_fn, func: apply_fn(lambda t: func(*t)),
        features=features,
        cpus=cpus,
    )
