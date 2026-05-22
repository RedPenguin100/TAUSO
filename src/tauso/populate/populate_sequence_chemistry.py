import logging
import time
from typing import Iterable, Optional, Tuple

import pandas as pd

logger = logging.getLogger(__name__)

from ..data.consts import CHEMICAL_PATTERN, SEQUENCE
from ..features.sequence.seq_chemistry import gap_gc_content, wing_gap_gc_delta
from ..parallel_utils import make_apply_fn
from ..timer import Timer

FEATURE_SPECS: list[tuple[str, callable]] = [
    ("SeqChem_gap_gc_content", gap_gc_content),
    ("SeqChem_wing_gap_gc_delta", wing_gap_gc_delta),
]


def calc_feature(df: pd.DataFrame, col_name: str, func, cpus: int = 1, verbose=False) -> None:
    """
    Computes a feature using both sequence and chemical pattern.
    Uses parallel_apply if available and cpus > 1.
    """
    start_time = time.time()

    logger.debug("Starting %s...", col_name)

    # use_memory_fs=False prevents silent hangs in containerized/Docker environments
    apply_fn = make_apply_fn(df, n_jobs=cpus, progress_bar=verbose, verbose=0, use_memory_fs=False)
    df[col_name] = apply_fn(lambda row: func(row[SEQUENCE], row[CHEMICAL_PATTERN]), axis=1)

    duration = time.time() - start_time
    logger.debug("Finished %s | Time: %.2fs", col_name, duration)


def populate_sequence_chemistry_features(
    df,
    features: Optional[Iterable[str]] = None,
    cpus: int = 1,
) -> Tuple:
    available = {name: fn for name, fn in FEATURE_SPECS}
    feature_names = list(features) if features is not None else [name for name, _ in FEATURE_SPECS]

    for name in feature_names:
        # Wrap each feature calculation in the Timer context manager
        with Timer(name):
            calc_feature(df, name, available[name], cpus=cpus)

    return df, feature_names
