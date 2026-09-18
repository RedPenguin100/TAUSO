"""An unrecognised cell line has no expression, and that is a NaN feature, not a crash."""

import numpy as np
import pandas as pd
import pytest

from tauso.data.consts import CANONICAL_GENE_NAME, CELL_LINE_DEPMAP
from tauso.populate.populate_context import (
    populate_special_gene_expression,
    populate_special_transcript_expression,
    populate_target_dominant_transcript,
    populate_target_expression,
)

ROWS = pd.DataFrame(
    {
        CELL_LINE_DEPMAP: ["ACH-999999", None],
        CANONICAL_GENE_NAME: ["MALAT1", "TP53"],
    }
)

POPULATE_FUNCTIONS = [
    populate_target_expression,
    populate_special_gene_expression,
    populate_target_dominant_transcript,
    populate_special_transcript_expression,
]


@pytest.mark.parametrize("populate", POPULATE_FUNCTIONS, ids=lambda f: f.__name__)
def test_no_expression_at_all_scores_nan(populate):
    """Nothing resolved: a missing pair is a missing pair, whether it is one of them or all."""
    data, features = populate(ROWS.copy(), {})

    assert len(data) == len(ROWS)
    assert np.isnan(data[features].to_numpy(dtype=float)).all()


@pytest.mark.parametrize("populate", POPULATE_FUNCTIONS, ids=lambda f: f.__name__)
def test_the_features_are_named_either_way(populate):
    """The columns exist whether or not anything filled them, so the step's contract holds."""
    data, features = populate(ROWS.copy(), {})

    assert features
    assert set(features) <= set(data.columns)
