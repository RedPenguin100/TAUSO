"""Shared driver for the ``(name, func)`` feature-spec pattern used by the populate modules."""

import logging
from collections.abc import Callable, Iterable, Sequence

import pandas as pd

from ..pandas_utils import add_columns
from ..parallel_utils import make_apply_fn
from ..timer import Timer

logger = logging.getLogger(__name__)

FeatureSpec = tuple[str, Callable]


def compute_features(
    df: pd.DataFrame,
    specs: Sequence[FeatureSpec],
    apply_target: "pd.DataFrame | pd.Series",
    apply_feature: Callable[[Callable, Callable], "pd.Series"],
    features: Iterable[str] | None = None,
    cpus: int = 1,
    verbose: bool = False,
) -> tuple[pd.DataFrame, list[str]]:
    """Compute each ``(name, func)`` in ``specs`` into ``df[name]``, logging per-feature timing.

    ``apply_target`` is what gets applied over in parallel: the DataFrame for row-based
    features, or a prepared Series for element-based ones. ``apply_feature(apply_fn, func)``
    returns the resulting column and owns how ``func`` reads its inputs — which row columns,
    or a bare element.
    """
    available = dict(specs)
    names = list(features) if features is not None else [name for name, _ in specs]
    apply_fn = make_apply_fn(apply_target, n_jobs=cpus, progress_bar=verbose, verbose=0, use_memory_fs=False)

    # Attached in one pass: a column at a time fragments the frame and pandas warns.
    with Timer(log=False) as timer:
        computed = {name: apply_feature(apply_fn, available[name]) for name in names}
    logger.info("Computed %d features in %.4fs: %s", len(names), timer.elapsed_time, ", ".join(names))

    return add_columns(df, computed), names
