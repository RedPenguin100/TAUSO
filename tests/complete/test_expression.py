from pathlib import Path

import pandas as pd
import pytest

from tauso.data.consts import CELL_LINE_DEPMAP
from tauso.populate.populate_context import (
    populate_special_gene_expression,
    populate_target_expression,
)

_FIXTURE = Path(__file__).parent / "test_expression" / "expression_fixture.parquet"


@pytest.fixture(scope="session")
def expression_transcriptomes():
    """
    Committed fixture of expression values for the test cell lines and genes.
    Decouples regression tests from live DepMap data so they are deterministic
    regardless of DepMap cache state or data updates.
    """
    df = pd.read_parquet(_FIXTURE)
    transcriptomes = {}
    for ach_id, group in df.groupby("ModelID"):
        transcriptomes[ach_id] = group[["Gene", "expression_norm"]].reset_index(drop=True)
    return transcriptomes


@pytest.mark.parametrize("mini_sampled_data", [1000], indirect=True)
def test_target_expression_regression(mini_sampled_data, expression_transcriptomes, dataframe_regression):
    data = mini_sampled_data.copy()
    result, feats = populate_target_expression(data, expression_transcriptomes)
    dataframe_regression.check(result[feats])


@pytest.mark.parametrize("mini_sampled_data", [1000], indirect=True)
def test_special_gene_expression_regression(mini_sampled_data, expression_transcriptomes, dataframe_regression):
    data = mini_sampled_data.copy()
    result, feats = populate_special_gene_expression(data, expression_transcriptomes)
    dataframe_regression.check(result[feats])
