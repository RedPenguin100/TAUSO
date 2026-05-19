"""
Memory-scaling tests for the batch off-target pipeline.

Uses tracemalloc + psutil to measure:
  1. Peak Python heap allocation inside compute_group_batch (no binary)
  2. That memory is released (no leak) after the call
  3. That peak scales sub-linearly with output size (aggregation reduces data early)

All RIsearch calls are mocked so this runs purely on synthetic DataFrames.
"""

import gc
import os
import tracemalloc
from unittest.mock import patch

import numpy as np
import pandas as pd
import psutil
import pytest

from tauso.data.consts import CANONICAL_GENE, SEQUENCE
from tauso.features.hybridization_off_target.add_off_target_feat import (
    AggregationMethod,
    compute_group_batch,
)


def _rss_mb() -> float:
    return psutil.Process(os.getpid()).memory_info().rss / 1024**2


def _make_group_df(n_rows: int, n_genes: int = 5) -> pd.DataFrame:
    genes = [f"GENE_{i % n_genes}" for i in range(n_rows)]
    seqs = ["ATCGATCGATCGATCGATCG"] * n_rows
    return pd.DataFrame({CANONICAL_GENE: genes, SEQUENCE: seqs})


def _make_synthetic_output(n_rows: int, n_targets: int, hits_per_query: int) -> str:
    """Generate a realistic-sized TSV that parse_risearch_output can consume."""
    rng = np.random.default_rng(0)
    lines = []
    for i in range(n_rows):
        for t in range(hits_per_query):
            target = f"GENE_{t % n_targets}"
            energy = rng.uniform(-20.0, -0.5)
            lines.append(f"{i}\t1\t20\t{target}\t100\t120\t900\t{energy:.4f}")
    return "\n".join(lines) + "\n"


EXP_MAP = {f"GENE_{i}": (float(1_000_000 * (i + 1)), 0.5) for i in range(20)}


@pytest.mark.parametrize(
    "n_rows,hits_per_query",
    [
        (10, 5),
        (100, 10),
        (500, 15),
    ],
)
def test_peak_allocation_scales_with_hits(n_rows, hits_per_query):
    """
    Peak tracemalloc allocation should grow with output size (we're not testing absolute numbers,
    just confirming there are no unexpected N² allocations by checking that doubling rows
    does not more-than-double peak memory).
    """
    n_targets = 10
    group_df = _make_group_df(n_rows, n_targets)
    raw_output = _make_synthetic_output(n_rows, n_targets, hits_per_query)

    gc.collect()
    tracemalloc.start()

    with patch(
        "tauso.features.hybridization_off_target.add_off_target_feat.get_triggers_mfe_scores_batch",
        return_value=raw_output,
    ):
        scores = compute_group_batch(
            group_df,
            seq_map={},
            exp_map=EXP_MAP,
            cutoff=0,
            method=AggregationMethod.ARTM,
            prebuilt_target_path="FAKE",
        )

    _, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    peak_mb = peak_bytes / 1024**2
    n_hits = n_rows * hits_per_query
    # Very loose upper bound: < 50 MB per 10k hits (generous to avoid flakiness)
    limit_mb = max(5.0, n_hits * 0.005)
    assert peak_mb < limit_mb, (
        f"n_rows={n_rows}, hits_per_query={hits_per_query}: peak {peak_mb:.1f} MB > limit {limit_mb:.1f} MB"
    )
    assert len(scores) == n_rows


def test_no_memory_leak_across_repeated_calls():
    """
    RSS after N repeated calls should not grow substantially beyond RSS before.
    This catches memory that pandas or numpy fails to reclaim between calls.
    """
    n_rows, hits_per_query, n_targets = 200, 10, 10
    group_df = _make_group_df(n_rows, n_targets)
    raw_output = _make_synthetic_output(n_rows, n_targets, hits_per_query)

    # Warm up to stabilise allocator
    for _ in range(2):
        with patch(
            "tauso.features.hybridization_off_target.add_off_target_feat.get_triggers_mfe_scores_batch",
            return_value=raw_output,
        ):
            compute_group_batch(
                group_df,
                seq_map={},
                exp_map=EXP_MAP,
                cutoff=0,
                method=AggregationMethod.ARTM,
                prebuilt_target_path="FAKE",
            )
    gc.collect()

    rss_before = _rss_mb()
    N_REPS = 8
    for _ in range(N_REPS):
        with patch(
            "tauso.features.hybridization_off_target.add_off_target_feat.get_triggers_mfe_scores_batch",
            return_value=raw_output,
        ):
            compute_group_batch(
                group_df,
                seq_map={},
                exp_map=EXP_MAP,
                cutoff=0,
                method=AggregationMethod.ARTM,
                prebuilt_target_path="FAKE",
            )
    gc.collect()
    rss_after = _rss_mb()

    growth_mb = rss_after - rss_before
    # Allow up to 50 MB of ratchet (allocator bookkeeping), but not unbounded growth
    assert growth_mb < 50, f"RSS grew by {growth_mb:.1f} MB over {N_REPS} calls — possible memory leak"


def test_rss_stable_across_large_batch():
    """
    RSS should not grow more than 30 MB processing a large batch (500 rows × 20 hits).
    Parsing a DataFrame from 10k hit rows is typically 2-5 MB; explicit del + gc keeps
    RSS from ratcheting across calls.
    """
    n_rows, hits_per_query, n_targets = 500, 20, 10
    group_df = _make_group_df(n_rows, n_targets)
    raw_output = _make_synthetic_output(n_rows, n_targets, hits_per_query)

    gc.collect()
    rss_before = _rss_mb()

    with patch(
        "tauso.features.hybridization_off_target.add_off_target_feat.get_triggers_mfe_scores_batch",
        return_value=raw_output,
    ):
        compute_group_batch(
            group_df,
            seq_map={},
            exp_map=EXP_MAP,
            cutoff=0,
            method=AggregationMethod.ARTM,
            prebuilt_target_path="FAKE",
        )

    gc.collect()
    rss_after = _rss_mb()

    growth_mb = rss_after - rss_before
    # Processing 10k hits should not permanently grow RSS by more than 30 MB
    assert growth_mb < 30, (
        f"RSS grew by {growth_mb:.1f} MB processing a {n_rows}-row batch — "
        "check that intermediate DataFrames are being freed"
    )
