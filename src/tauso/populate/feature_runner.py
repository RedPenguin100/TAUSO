"""Shared driver for the ``(name, func)`` feature-spec pattern used by the populate modules."""

from collections.abc import Callable, Iterable, Sequence

import pandas as pd

from ..parallel_utils import make_apply_fn
from ..timer import Timer

FeatureSpec = tuple[str, Callable]


def run_feature_specs(
    df: pd.DataFrame,
    specs: Sequence[FeatureSpec],
    apply_target: "pd.DataFrame | pd.Series",
    apply_feature: Callable[[Callable, Callable], "pd.Series"],
    features: Iterable[str] | None = None,
    cpus: int = 1,
    verbose: bool = False,
) -> tuple[pd.DataFrame, list[str]]:
    """Compute each ``(name, func)`` in ``specs`` into ``df[name]``, timed per feature.

    ``apply_target`` is what gets applied over in parallel: the DataFrame for row-based
    features, or a prepared Series for element-based ones. ``apply_feature(apply_fn, func)``
    returns the resulting column and owns how ``func`` reads its inputs — which row columns,
    or a bare element.
    """
    available = dict(specs)
    names = list(features) if features is not None else [name for name, _ in specs]
    apply_fn = make_apply_fn(apply_target, n_jobs=cpus, progress_bar=verbose, verbose=0, use_memory_fs=False)
    for name in names:
        with Timer(name):
            df[name] = apply_feature(apply_fn, available[name])
    return df, names
