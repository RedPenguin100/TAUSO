import logging
from typing import Iterable, Optional, Tuple

logger = logging.getLogger(__name__)

from ..data.consts import ASO_SEQUENCE, CHEMICAL_PATTERN
from ..features.sequence.seq_chemistry import gap_gc_content, wing_gap_gc_delta
from .feature_runner import run_feature_specs

FEATURE_SPECS: list[tuple[str, callable]] = [
    ("mod_sugar_gap_gc_content", gap_gc_content),
    ("mod_sugar_wing_gap_gc_delta", wing_gap_gc_delta),
]


def populate_sequence_chemistry_features(
    df,
    features: Optional[Iterable[str]] = None,
    cpus: int = 1,
) -> Tuple:
    return run_feature_specs(
        df,
        FEATURE_SPECS,
        df,
        lambda apply_fn, func: apply_fn(lambda row: func(row[ASO_SEQUENCE], row[CHEMICAL_PATTERN]), axis=1),
        features=features,
        cpus=cpus,
    )
