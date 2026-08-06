import pytest

from notebooks.data.OligoAI.parse_chemistry import assign_chemistry
from tauso.populate.populate_modification import populate_modifications
from tests.complete.conftest import get_n_jobs


@pytest.fixture(scope="session")
def chemistry_only_data(base_data):
    """Modification features only consume CHEMICAL_PATTERN. assign_chemistry derives it
    directly from MODIFICATION_STRING, so this fixture skips the structure / genome-DB stack
    that other tests pull in via `chemistry_data`.
    """
    return assign_chemistry(base_data)


@pytest.fixture(scope="session")
def mini_sampled_chemistry_only_data(request, chemistry_only_data):
    n_samples = getattr(request, "param", 1000)
    actual_samples = min(n_samples, len(chemistry_only_data))
    return chemistry_only_data.sample(n=actual_samples, random_state=42).copy()


@pytest.mark.parametrize("mini_sampled_chemistry_only_data", [10000], indirect=True)
def test_modification_sugar(mini_sampled_chemistry_only_data, dataframe_regression):
    data, feature_names = populate_modifications(mini_sampled_chemistry_only_data, n_cores=get_n_jobs())
    dataframe_regression.check(data[["index_oligo"] + feature_names])
