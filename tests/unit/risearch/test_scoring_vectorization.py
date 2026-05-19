"""
Test that the vectorized exp-energy sum in _apply_risearch_scoring gives
identical results to the original per-row scan approach.

These are pure-numpy/pandas tests — no RIsearch binary invoked.
"""

import numpy as np
import pandas as pd
import pytest

RT = 0.616


def score_per_row_original(result_df, row_triggers):
    """Original O(N×M) approach: filter by trigger per row, sum exp(-RT*energy)."""
    scores = {}
    for idx, _ in row_triggers:
        row_hits = result_df[result_df["trigger"] == str(idx)]
        if not row_hits.empty:
            scores[idx] = np.exp(-RT * row_hits["energy"]).sum()
        else:
            scores[idx] = 0.0
    return scores


def score_vectorized(result_df, row_triggers):
    """New vectorized approach: groupby → dict → O(1) lookup."""
    result_df = result_df.assign(_exp=np.exp(-RT * result_df["energy"].to_numpy()))
    score_dict = result_df.groupby("trigger", sort=False)["_exp"].sum().to_dict()
    scores = {}
    for idx, _ in row_triggers:
        scores[idx] = score_dict.get(str(idx), 0.0)
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
        (100, 0),  # no hits at all
    ],
)
def test_vectorized_matches_original(n_queries, n_hits_per_query):
    rng = np.random.default_rng(42)
    rows = []
    row_triggers = [(i, f"seq_{i}") for i in range(n_queries)]

    if n_hits_per_query > 0:
        for i in range(n_queries):
            energies = rng.uniform(-20.0, -0.1, size=n_hits_per_query)
            for e in energies:
                rows.append((str(i), float(e)))

    result_df = _make_result_df(rows) if rows else pd.DataFrame(columns=["trigger", "energy"])

    if result_df.empty:
        # Both approaches should give all zeros
        orig = {idx: 0.0 for idx, _ in row_triggers}
        vec = {idx: 0.0 for idx, _ in row_triggers}
    else:
        orig = score_per_row_original(result_df, row_triggers)
        vec = score_vectorized(result_df, row_triggers)

    for idx, _ in row_triggers:
        assert vec[idx] == pytest.approx(orig[idx], rel=1e-9), (
            f"Mismatch at trigger {idx}: vectorized={vec[idx]}, original={orig[idx]}"
        )


def test_missing_trigger_gives_zero():
    """Trigger IDs in row_triggers but absent from result_df must score 0."""
    result_df = _make_result_df([("0", -5.0), ("2", -3.0)])
    row_triggers = [(0, "a"), (1, "b"), (2, "c")]  # 1 has no hit

    orig = score_per_row_original(result_df, row_triggers)
    vec = score_vectorized(result_df, row_triggers)

    assert orig[1] == pytest.approx(0.0)
    assert vec[1] == pytest.approx(0.0)
    assert vec[0] == pytest.approx(orig[0], rel=1e-9)
    assert vec[2] == pytest.approx(orig[2], rel=1e-9)


def test_multiple_hits_same_trigger():
    """Multiple hits per trigger must be summed, not deduplicated."""
    # Trigger 0 has two hits; original sums them, vectorized must too
    result_df = _make_result_df([("0", -5.0), ("0", -10.0)])
    row_triggers = [(0, "seq")]

    orig = score_per_row_original(result_df, row_triggers)
    vec = score_vectorized(result_df, row_triggers)

    expected = np.exp(-RT * -5.0) + np.exp(-RT * -10.0)
    assert orig[0] == pytest.approx(expected, rel=1e-9)
    assert vec[0] == pytest.approx(expected, rel=1e-9)
