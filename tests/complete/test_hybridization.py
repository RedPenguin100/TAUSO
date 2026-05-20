import pytest

from tauso.populate.populate_hybridization import populate_hybridization


@pytest.fixture(scope="session")
def mini_sampled_data(request, chemistry_data):
    """Samples chemistry_data once per (n_samples, session).

    chemistry_data is the session-scoped fixture from conftest (built on structure_data).
    Tests must call .copy() before mutating the returned DataFrame.
    """
    n_samples = getattr(request, "param", 1000)
    actual_samples = min(n_samples, len(chemistry_data))
    return chemistry_data.sample(n=actual_samples, random_state=42).copy()


@pytest.mark.parametrize("mini_sampled_data", [10000], indirect=True)
def test_hybridization(mini_sampled_data, dataframe_regression):
    # populate_hybridization copies internally — no need to copy here
    data, feature_names = populate_hybridization(mini_sampled_data)
    dataframe_regression.check(data[["index_oligo"] + feature_names])
