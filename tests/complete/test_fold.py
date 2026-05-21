import gc
import time
import tracemalloc

import pytest
import logging

from tauso.populate.populate_fold import (
    populate_mfe_features,
    populate_sense_accessibility_batch,
)


ACCESS_CONFIGURATIONS = [
    {"flank": 120, "access": 20, "seeds": [13]},
    {"flank": 120, "access": 13, "seeds": [4, 6, 8]},
    {"flank": 120, "access": 20, "seeds": [4, 6, 8]},
    {"flank": 120, "access": 13, "seeds": [13, 26, 39]},
]


@pytest.mark.parametrize("mini_sampled_data", [1000], indirect=True)
@pytest.mark.parametrize("config", ACCESS_CONFIGURATIONS)
def test_access(mini_sampled_data, gene_to_data, config, dataframe_regression):
    data, feature_name = populate_sense_accessibility_batch(
        mini_sampled_data,
        gene_to_data,
        batch_size=250,
        flank_size=config["flank"],
        access_size=config["access"],
        seed_sizes=config["seeds"],
        n_jobs=4,
    )
    dataframe_regression.check(data[["index_oligo", feature_name]])


@pytest.mark.parametrize("mini_sampled_data", [1000], indirect=True)
def test_mfe(mini_sampled_data, gene_to_data, dataframe_regression):
    # populate_mfe_features mutates its input — copy to protect the session fixture
    data, feature_names = populate_mfe_features(mini_sampled_data.copy(), gene_to_data)
    dataframe_regression.check(data[["index_oligo"] + feature_names])


@pytest.mark.parametrize("mini_sampled_data", [1000], indirect=True)
def test_mfe_parallel_matches_baseline(mini_sampled_data, gene_to_data, dataframe_regression):
    """n_jobs=4 output must be numerically identical to the single-worker baseline."""
    data, feature_names = populate_mfe_features(mini_sampled_data.copy(), gene_to_data, n_jobs=4)
    # Reuse the existing n_jobs=1 baseline — same results expected
    dataframe_regression.check(data[["index_oligo"] + feature_names], basename="test_mfe_1000_")


@pytest.mark.parametrize("mini_sampled_data", [1000], indirect=True)
def test_mfe_parallel_benchmark(mini_sampled_data, gene_to_data):
    """Print wall time and peak Python heap for n_jobs=1 vs n_jobs=4."""
    df = mini_sampled_data.copy()

    def _run(n_jobs):
        gc.collect()
        tracemalloc.start()
        t0 = time.perf_counter()
        populate_mfe_features(df.copy(), gene_to_data, n_jobs=n_jobs)
        elapsed = time.perf_counter() - t0
        _, peak_bytes = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        return elapsed, peak_bytes / 1024**2

    t1, mem1 = _run(1)
    t4, mem4 = _run(4)

    print(f"\n{'':>10}  {'time':>8}  {'peak heap':>10}")
    print(f"{'n_jobs=1':>10}  {t1:>7.1f}s  {mem1:>9.1f}MB")
    print(f"{'n_jobs=4':>10}  {t4:>7.1f}s  {mem4:>9.1f}MB")
    print(f"{'speedup':>10}  {t1/t4:>8.2f}x")
