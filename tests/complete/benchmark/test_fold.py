import gc
import time
import tracemalloc

import pytest

from tauso.populate.populate_fold import populate_mfe_features


@pytest.mark.benchmark
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
