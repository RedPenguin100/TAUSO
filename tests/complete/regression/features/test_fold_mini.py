"""Mini MFE fold regression test at 100 samples.

Wired to the module-level defaults of populate_fold so that re-parametrizing
the defaults forces a baseline regen. Also reports wall-clock so we can
compare runtime across parameter sets.
"""

import time

import pytest
from tests.complete.conftest import get_n_jobs

from tauso.populate.populate_fold import populate_mfe_features


@pytest.mark.parametrize("mini_sampled_data", [100], indirect=True)
def test_mfe_mini(mini_sampled_data, gene_to_data, dataframe_regression):
    t0 = time.perf_counter()
    data, feature_names = populate_mfe_features(
        mini_sampled_data.copy(),
        gene_to_data,
        n_jobs=get_n_jobs(),
    )
    elapsed = time.perf_counter() - t0
    print(
        f"\n[MFE mini] n_jobs={get_n_jobs()} "
        f"n_rows={len(mini_sampled_data)} "
        f"features={len(feature_names)} "
        f"elapsed={elapsed:.2f}s "
        f"per_row={elapsed / len(mini_sampled_data) * 1000:.1f}ms"
    )
    dataframe_regression.check(data[["index_oligo"] + feature_names])
