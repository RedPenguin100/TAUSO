"""
Verify that _sum_exp_energy_by_trigger (the production vectorized function)
gives results identical to a naive per-row reference implementation.
"""

import numpy as np
import pandas as pd
import pytest

from tauso.features.hybridization.off_target.off_target_specific_gene import (
    _RT,
    _sum_exp_energy_by_trigger,
)


def _score_per_row_reference(result_df, row_triggers):
    """Naive O(N×M) reference: filter hits per trigger, sum exp(-RT*energy)."""
    scores = {}
    for idx, _ in row_triggers:
        row_hits = result_df[result_df["trigger"] == str(idx)]
        scores[idx] = float(np.exp(-_RT * row_hits["energy"]).sum()) if not row_hits.empty else 0.0
    return scores


def _make_result_df(rows):
    """rows: list of (trigger_id_str, energy)"""
    return pd.DataFrame(rows, columns=["trigger", "energy"])


@pytest.mark.parametrize(
    "n_queries,n_hits_per_query",
    [
        (1, 1),
        (5, 3),
        (50, 10),
        (100, 0),
    ],
)
def test_production_matches_reference(n_queries, n_hits_per_query):
    rng = np.random.default_rng(42)
    rows = []
    row_triggers = [(i, f"seq_{i}") for i in range(n_queries)]

    if n_hits_per_query > 0:
        for i in range(n_queries):
            for e in rng.uniform(-20.0, -0.1, size=n_hits_per_query):
                rows.append((str(i), float(e)))

    if not rows:
        return  # nothing to compare for zero-hit case

    result_df = _make_result_df(rows)
    ref = _score_per_row_reference(result_df, row_triggers)
    score_dict = _sum_exp_energy_by_trigger(result_df)
    prod = {idx: score_dict.get(str(idx), 0.0) for idx, _ in row_triggers}

    for idx, _ in row_triggers:
        assert prod[idx] == pytest.approx(ref[idx], rel=1e-9), (
            f"Mismatch at trigger {idx}: production={prod[idx]}, reference={ref[idx]}"
        )


def test_missing_trigger_gives_zero():
    result_df = _make_result_df([("0", -5.0), ("2", -3.0)])
    row_triggers = [(0, "a"), (1, "b"), (2, "c")]

    ref = _score_per_row_reference(result_df, row_triggers)
    score_dict = _sum_exp_energy_by_trigger(result_df)
    prod = {idx: score_dict.get(str(idx), 0.0) for idx, _ in row_triggers}

    assert prod[1] == pytest.approx(0.0)
    assert ref[1] == pytest.approx(0.0)
    assert prod[0] == pytest.approx(ref[0], rel=1e-9)
    assert prod[2] == pytest.approx(ref[2], rel=1e-9)


def test_multiple_hits_same_trigger():
    result_df = _make_result_df([("0", -5.0), ("0", -10.0)])
    row_triggers = [(0, "seq")]

    ref = _score_per_row_reference(result_df, row_triggers)
    score_dict = _sum_exp_energy_by_trigger(result_df)
    prod = {idx: score_dict.get(str(idx), 0.0) for idx, _ in row_triggers}

    expected = np.exp(-_RT * -5.0) + np.exp(-_RT * -10.0)
    assert ref[0] == pytest.approx(expected, rel=1e-9)
    assert prod[0] == pytest.approx(expected, rel=1e-9)
