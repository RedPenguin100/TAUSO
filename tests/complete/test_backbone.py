import pandas as pd
import pytest

from notebooks.consts import OLIGO_CSV_INDEXED
from notebooks.data.OligoAI.parse_chemistry import assign_chemistry
from notebooks.data.OligoAI.utility import standardize_oligo_ai_data
from tauso.populate.populate_backbone import populate_backbone_features


@pytest.fixture(scope="module")
def chemistry_sample(request):
    n = getattr(request, "param", 10000)
    data = assign_chemistry(standardize_oligo_ai_data(pd.read_csv(OLIGO_CSV_INDEXED, engine="pyarrow")))
    return data.sample(n=min(n, len(data)), random_state=42).reset_index(drop=True)


@pytest.mark.parametrize("chemistry_sample", [10000], indirect=True)
def test_backbone_features(chemistry_sample, dataframe_regression):
    data = chemistry_sample.copy()
    processed_data, feature_cols = populate_backbone_features(data)
    dataframe_regression.check(processed_data[feature_cols])
