import pytest

from notebooks.data.OligoAI.parse_chemistry import (
    assign_chemistry,
)
from tauso.populate.populate_rnase import populate_rnase_features


@pytest.fixture(scope="module")
def chemistry_data(base_data):
    return assign_chemistry(base_data)


@pytest.fixture
def sampled_chem_data(request, chemistry_data):
    n_samples = getattr(request, "param", 1000)
    actual_samples = min(n_samples, len(chemistry_data))
    return chemistry_data.sample(n=actual_samples, random_state=42).copy()


@pytest.mark.parametrize("sampled_chem_data", [10000], indirect=True)
def test_rnase_features_regression(sampled_chem_data, dataframe_regression):
    data = sampled_chem_data.copy()

    processed_data, rnase_features = populate_rnase_features(data)

    # Assuming rnase_features is a list of the new column names
    dataframe_regression.check(processed_data[rnase_features])
