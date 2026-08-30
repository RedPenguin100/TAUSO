"""Rows with no DepMap cell line must score NaN, audibly.

`populate_off_target_specific` groups by CELL_LINE_DEPMAP. pandas drops NaN keys, so when every
row was unresolved the group list came out empty and `pd.concat` raised "No objects to
concatenate" — contradicting the function's own docstring, which promises unresolved rows score
NaN. The fillna sentinel fixed the crash; these tests pin that, and pin that the all-NaN result
announces itself instead of looking like a successful run.
"""

import logging

import numpy as np
import pandas as pd
import pytest

from tauso.data.consts import CELL_LINE_DEPMAP
from tauso.populate.populate_off_target import UNRESOLVED_CELL_LINE, populate_off_target_specific


def _frame(depmap_ids):
    return pd.DataFrame(
        {
            "index_oligo": np.arange(len(depmap_ids)),
            CELL_LINE_DEPMAP: depmap_ids,
            "aso_sequence": ["ACGTACGTACGTACGTACGT"] * len(depmap_ids),
        }
    )


def _run(df, caplog):
    """No reference data, so every cell line is unknown and every row takes the NaN branch."""
    with caplog.at_level(logging.WARNING):
        return populate_off_target_specific(
            ASO_df=df,
            gene_to_data={},
            cell_line2data={},
            top_n_list=[50],
            cutoff_list=[1200],
            method="BOLTZ",
            n_jobs=1,
        )


def test_every_row_unresolved_scores_nan_instead_of_raising(caplog):
    df = _frame([None, None, None])
    out, feature_names = _run(df, caplog)

    assert len(out) == 3, "rows must be kept, not dropped"
    for col in feature_names:
        assert out[col].isna().all()


def test_every_row_unresolved_says_so(caplog):
    _run(_frame([None, None, None]), caplog)
    assert "no row has a resolvable DepMap cell line" in caplog.text


def test_some_rows_unresolved_reports_the_count(caplog):
    _run(_frame(["ACH-000681", None, None]), caplog)
    assert "2 of 3 rows have no DepMap cell line" in caplog.text


def test_all_resolved_says_nothing(caplog):
    _run(_frame(["ACH-000681", "ACH-000681"]), caplog)
    assert "no DepMap cell line" not in caplog.text


@pytest.mark.parametrize("ids", [[None], ["ACH-000681"]])
def test_the_sentinel_never_leaks_into_the_output(ids, caplog):
    out, _ = _run(_frame(ids), caplog)
    assert UNRESOLVED_CELL_LINE not in out[CELL_LINE_DEPMAP].astype(str).tolist()
