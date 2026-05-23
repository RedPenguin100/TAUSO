"""
Benchmark tests for all three off-target RIsearch pathways.

Measures wall time and peak RSS across:
  - n_jobs ∈ {1, 4, 8}
  - chunk_size ∈ {50, 100, 250, 500}
  - executor ∈ {"thread", "process"}

Run with:
    pytest tests/complete/benchmark/test_off_target_parallel.py \
        --run-benchmarks -s --n-jobs=4

The --n-jobs flag controls the *default* jobs used in complete/ session fixtures;
each benchmark also sweeps its own values independently.
"""

import gc
import os
import time
import tracemalloc

import psutil
import pytest

from tauso.features.hybridization_off_target.add_off_target_feat import AggregationMethod
from tauso.features.hybridization_off_target.off_target_feature import (
    populate_off_target_general,
    populate_off_target_specific,
)
from tauso.features.hybridization_off_target.off_target_specific_gene import (
    off_target_specific_seq_pandarallel,
    on_target_total_hybridization,
)

CUTOFFS = [800]
TOP_N = [25]
METHOD = AggregationMethod.ARTM

# These genes are always available via gene_to_data fixture
SINGLE_TARGET_GENES = ["RNASEH1", "ACTB"]


def _rss_mb() -> float:
    return psutil.Process(os.getpid()).memory_info().rss / 1024**2


def _timed_run(fn, *args, **kwargs):
    """Run fn(*args, **kwargs), return (result, elapsed_s, peak_heap_MB, rss_delta_MB)."""
    gc.collect()
    rss_before = _rss_mb()
    tracemalloc.start()
    t0 = time.perf_counter()
    result = fn(*args, **kwargs)
    elapsed = time.perf_counter() - t0
    _, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    gc.collect()
    rss_delta = _rss_mb() - rss_before
    return result, elapsed, peak_bytes / 1024**2, rss_delta


def _print_table(title, rows):
    """rows: list of (label, time_s, heap_MB, rss_delta_MB)"""
    print(f"\n{'─'*60}")
    print(f"  {title}")
    print(f"{'─'*60}")
    print(f"  {'config':<25} {'time':>8}  {'heap':>9}  {'rss Δ':>9}")
    print(f"  {'-'*25} {'------':>8}  {'-------':>9}  {'-------':>9}")
    for label, t, heap, rss in rows:
        print(f"  {label:<25} {t:>7.1f}s  {heap:>8.1f}M  {rss:>8.1f}M")
    print(f"{'─'*60}")


# ─────────────────────────────────────────────────────────────
# 1. Single-gene pathway (off_target_specific_seq_pandarallel)
# ─────────────────────────────────────────────────────────────

@pytest.mark.benchmark
@pytest.mark.parametrize("mini_structure_data", [500], indirect=True)
def test_single_gene_njobs_sweep(mini_structure_data, gene_to_data_full):
    """Wall time and memory for off_target_specific_seq vs n_jobs={1,4,8}."""
    data = mini_structure_data.copy()
    gene = SINGLE_TARGET_GENES[0]
    cutoff = CUTOFFS[0]

    rows = []
    for nj in (1, 4, 8):
        _, t, heap, rss = _timed_run(
            off_target_specific_seq_pandarallel, data, gene, gene_to_data_full, cutoff, n_jobs=nj
        )
        rows.append((f"n_jobs={nj}", t, heap, rss))

    _print_table("off_target_specific_seq: n_jobs sweep", rows)
    base_time = rows[0][1]
    for label, t, _, _ in rows[1:]:
        assert t < base_time * 1.5 or t < 5.0, f"{label} slower than serial by >50%: {t:.1f}s vs {base_time:.1f}s"


@pytest.mark.benchmark
@pytest.mark.parametrize("mini_structure_data", [500], indirect=True)
def test_single_gene_chunk_size_sweep(mini_structure_data, gene_to_data_full):
    """Wall time and memory for off_target_specific_seq vs chunk_size={50,100,250,500}."""
    data = mini_structure_data.copy()
    gene = SINGLE_TARGET_GENES[0]
    cutoff = CUTOFFS[0]

    rows = []
    for cs in (50, 100, 250, 500):
        _, t, heap, rss = _timed_run(
            off_target_specific_seq_pandarallel, data, gene, gene_to_data_full, cutoff, n_jobs=4, chunk_size=cs
        )
        rows.append((f"chunk_size={cs}", t, heap, rss))

    _print_table("off_target_specific_seq: chunk_size sweep (n_jobs=4)", rows)


@pytest.mark.benchmark
@pytest.mark.parametrize("mini_structure_data", [500], indirect=True)
def test_on_target_njobs_sweep(mini_structure_data, gene_to_data):
    """Wall time and memory for on_target_total_hybridization vs n_jobs={1,4,8}."""
    data = mini_structure_data.copy()
    cutoff = CUTOFFS[0]

    rows = []
    for nj in (1, 4, 8):
        _, t, heap, rss = _timed_run(
            on_target_total_hybridization, data, gene_to_data, cutoff, n_jobs=nj
        )
        rows.append((f"n_jobs={nj}", t, heap, rss))

    _print_table("on_target_total_hybridization: n_jobs sweep", rows)


@pytest.mark.benchmark
@pytest.mark.parametrize("mini_structure_data", [500], indirect=True)
def test_on_target_chunk_size_sweep(mini_structure_data, gene_to_data):
    """Wall time and memory for on_target_total_hybridization vs chunk_size={50,100,250,500}."""
    data = mini_structure_data.copy()
    cutoff = CUTOFFS[0]

    rows = []
    for cs in (50, 100, 250, 500):
        _, t, heap, rss = _timed_run(
            on_target_total_hybridization, data, gene_to_data, cutoff, n_jobs=4, chunk_size=cs
        )
        rows.append((f"chunk_size={cs}", t, heap, rss))

    _print_table("on_target_total_hybridization: chunk_size sweep (n_jobs=4)", rows)


# ─────────────────────────────────────────────────────────────
# 2. populate_off_target_specific  (cell-line pathway)
# ─────────────────────────────────────────────────────────────

@pytest.mark.benchmark
@pytest.mark.parametrize("mini_structure_data", [500], indirect=True)
def test_ot_specific_njobs_sweep(mini_structure_data, gene_to_data, transcriptomes):
    """Wall time and memory for populate_off_target_specific vs n_jobs={1,4,8}."""
    data = mini_structure_data.copy()

    rows = []
    for nj in (1, 4, 8):
        _, t, heap, rss = _timed_run(
            populate_off_target_specific,
            data, gene_to_data, transcriptomes, TOP_N, CUTOFFS, METHOD,
            n_jobs=nj,
        )
        rows.append((f"n_jobs={nj}", t, heap, rss))

    _print_table("populate_off_target_specific: n_jobs sweep", rows)


@pytest.mark.benchmark
@pytest.mark.parametrize("mini_structure_data", [500], indirect=True)
def test_ot_specific_chunk_size_sweep(mini_structure_data, gene_to_data, transcriptomes):
    """Wall time and memory for populate_off_target_specific vs chunk_size={50,100,250,500}."""
    data = mini_structure_data.copy()

    rows = []
    for cs in (50, 100, 250, 500):
        _, t, heap, rss = _timed_run(
            populate_off_target_specific,
            data, gene_to_data, transcriptomes, TOP_N, CUTOFFS, METHOD,
            n_jobs=4, chunk_size=cs,
        )
        rows.append((f"chunk_size={cs}", t, heap, rss))

    _print_table("populate_off_target_specific: chunk_size sweep (n_jobs=4)", rows)


@pytest.mark.benchmark
@pytest.mark.parametrize("mini_structure_data", [500], indirect=True)
def test_ot_specific_executor_sweep(mini_structure_data, gene_to_data, transcriptomes):
    """Compare ThreadPoolExecutor vs ProcessPoolExecutor for populate_off_target_specific."""
    data = mini_structure_data.copy()

    rows = []
    for exc in ("thread", "process"):
        _, t, heap, rss = _timed_run(
            populate_off_target_specific,
            data, gene_to_data, transcriptomes, TOP_N, CUTOFFS, METHOD,
            n_jobs=4, chunk_size=250, executor=exc,
        )
        rows.append((f"executor={exc}", t, heap, rss))

    _print_table("populate_off_target_specific: executor sweep (n_jobs=4, chunk=250)", rows)


# ─────────────────────────────────────────────────────────────
# 3. populate_off_target_general  (general pathway)
# ─────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def transcriptomes_with_general_bm(request):
    """Module-scoped version of the general transcriptome fixture for benchmarks."""
    from pathlib import Path

    from tauso.common.gtf import filter_gtf_genes
    from tauso.data.data import get_data_dir, load_gtf_db
    from tauso.features.hybridization_off_target.common import get_general_expression_of_genes

    transcriptomes = request.getfixturevalue("transcriptomes")
    db = load_gtf_db()
    valid_genes = filter_gtf_genes(db, filter_mode="non_mt")

    data_dir = get_data_dir()
    exp_path = Path(data_dir) / "OmicsExpressionTPMLogp1HumanAllGenesStranded.csv"

    mean_exp_data = get_general_expression_of_genes(exp_path, valid_genes)
    result = dict(transcriptomes)
    result["general"] = mean_exp_data
    return result


@pytest.mark.benchmark
@pytest.mark.parametrize("mini_structure_data", [500], indirect=True)
def test_ot_general_njobs_sweep(mini_structure_data, gene_to_data_full, transcriptomes_with_general_bm):
    """Wall time and memory for populate_off_target_general vs n_jobs={1,4,8}."""
    data = mini_structure_data.copy()

    rows = []
    for nj in (1, 4, 8):
        _, t, heap, rss = _timed_run(
            populate_off_target_general,
            data, gene_to_data_full, transcriptomes_with_general_bm, TOP_N, CUTOFFS, METHOD,
            n_jobs=nj,
        )
        rows.append((f"n_jobs={nj}", t, heap, rss))

    _print_table("populate_off_target_general: n_jobs sweep", rows)


@pytest.mark.benchmark
@pytest.mark.parametrize("mini_structure_data", [500], indirect=True)
def test_ot_general_chunk_size_sweep(mini_structure_data, gene_to_data_full, transcriptomes_with_general_bm):
    """Wall time and memory for populate_off_target_general vs chunk_size={50,100,250,500}."""
    data = mini_structure_data.copy()

    rows = []
    for cs in (50, 100, 250, 500):
        _, t, heap, rss = _timed_run(
            populate_off_target_general,
            data, gene_to_data_full, transcriptomes_with_general_bm, TOP_N, CUTOFFS, METHOD,
            n_jobs=4, chunk_size=cs,
        )
        rows.append((f"chunk_size={cs}", t, heap, rss))

    _print_table("populate_off_target_general: chunk_size sweep (n_jobs=4)", rows)


@pytest.mark.benchmark
@pytest.mark.parametrize("mini_structure_data", [500], indirect=True)
def test_ot_general_executor_sweep(mini_structure_data, gene_to_data_full, transcriptomes_with_general_bm):
    """Compare ThreadPoolExecutor vs ProcessPoolExecutor for populate_off_target_general."""
    data = mini_structure_data.copy()

    rows = []
    for exc in ("thread", "process"):
        _, t, heap, rss = _timed_run(
            populate_off_target_general,
            data, gene_to_data_full, transcriptomes_with_general_bm, TOP_N, CUTOFFS, METHOD,
            n_jobs=4, chunk_size=250, executor=exc,
        )
        rows.append((f"executor={exc}", t, heap, rss))

    _print_table("populate_off_target_general: executor sweep (n_jobs=4, chunk=250)", rows)


# ─────────────────────────────────────────────────────────────
# 4. Combined sweep: best config vs worst
# ─────────────────────────────────────────────────────────────

@pytest.mark.benchmark
@pytest.mark.parametrize("mini_structure_data", [1000], indirect=True)
def test_ot_specific_combined_best_vs_worst(mini_structure_data, gene_to_data, transcriptomes):
    """
    Direct comparison of n_jobs=1 (worst) vs the empirically best (n_jobs=8, chunk=250).
    Uses 1000 ASOs to amplify differences.
    """
    data = mini_structure_data.copy()

    configs = [
        ("n_jobs=1 (serial)",     dict(n_jobs=1, chunk_size=250)),
        ("n_jobs=4 chunk=100",    dict(n_jobs=4, chunk_size=100)),
        ("n_jobs=4 chunk=250",    dict(n_jobs=4, chunk_size=250)),
        ("n_jobs=8 chunk=250",    dict(n_jobs=8, chunk_size=250)),
        ("n_jobs=8 chunk=500",    dict(n_jobs=8, chunk_size=500)),
    ]

    rows = []
    for label, kw in configs:
        _, t, heap, rss = _timed_run(
            populate_off_target_specific,
            data, gene_to_data, transcriptomes, TOP_N, CUTOFFS, METHOD,
            **kw,
        )
        rows.append((label, t, heap, rss))

    _print_table("populate_off_target_specific: combined sweep (1000 ASOs)", rows)
    serial_time = rows[0][1]
    print(f"\n  Speedups vs n_jobs=1:")
    for label, t, _, _ in rows[1:]:
        print(f"    {label:<30} {serial_time / t:>5.2f}×")
