import pytest
from tauso.populate.populate_hybridization import populate_hybridization



@pytest.fixture(scope="session")
def mini_chemistry_data(request, chemistry_data):
    """Samples the fully processed DataFrame right before the test runs."""
    n_samples = getattr(request, "param", 1000)
    actual_samples = min(n_samples, len(chemistry_data))
    return chemistry_data.sample(n=actual_samples, random_state=42).copy()


@pytest.mark.parametrize("mini_chemistry_data", [10000], indirect=True)
def test_hybridization(mini_chemistry_data, dataframe_regression):
    mini_chemistry_data, feature_names = populate_hybridization(mini_chemistry_data)

    dataframe_regression.check(mini_chemistry_data[["index_oligo"] + feature_names])
