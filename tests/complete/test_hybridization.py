import pytest

from notebooks.data.OligoAI.parse_chemistry import populate_chemistry
from tauso.populate.populate_hybridization import populate_hybridization


@pytest.fixture(scope="module")
def chemistry_data(final_data):
    return populate_chemistry(final_data)


@pytest.fixture
def mini_sampled_data(request, chemistry_data):
    """Samples the fully processed DataFrame right before the test runs."""
    n_samples = getattr(request, "param", 1000)
    actual_samples = min(n_samples, len(chemistry_data))
    return chemistry_data.sample(n=actual_samples, random_state=42).copy()


@pytest.mark.parametrize("mini_sampled_data", [10000], indirect=True)
def test_hybridization(mini_sampled_data, dataframe_regression):
    mini_sampled_data, feature_names = populate_hybridization(mini_sampled_data)

    dataframe_regression.check(mini_sampled_data[["index_oligo"] + feature_names])
