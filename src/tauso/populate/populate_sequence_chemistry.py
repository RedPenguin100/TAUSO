import time
import pandas as pd

from typing import Iterable, Optional, Tuple

from ..data.consts import SEQUENCE, CHEMICAL_PATTERN
from ..features.sequence.seq_chemistry import gap_gc_content, wing_gap_gc_delta
from ..timer import Timer

FEATURE_SPECS: list[tuple[str, callable]] = [
    ("SeqChem_gap_gc_content", gap_gc_content),
    ("SeqChem_wing_gap_gc_delta", wing_gap_gc_delta),
]


def calc_feature(
        df: pd.DataFrame,
        col_name: str,
        func,
        cpus: int = 1, verbose=False
) -> None:
    """
    Computes a feature using both sequence and chemical pattern.
    Uses parallel_apply if available.
    """
    start_time = time.time()

    if verbose:
        print(f"► Starting {col_name}...")

    try:
        from pandarallel import pandarallel
        pandarallel.initialize(progress_bar=verbose, verbose=0, nb_workers=cpus)

        # Apply across rows (axis=1) to pass both columns to the function
        df[col_name] = df.parallel_apply(
            lambda row: func(row[SEQUENCE], row[CHEMICAL_PATTERN]),
            axis=1
        )
    except Exception:
        # Fallback to standard pandas apply
        df[col_name] = df.apply(
            lambda row: func(row[SEQUENCE], row[CHEMICAL_PATTERN]),
            axis=1
        )

    if verbose:
        duration = time.time() - start_time
        print(f"✔ Finished {col_name} | Time: {duration:.2f}s\n")

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
