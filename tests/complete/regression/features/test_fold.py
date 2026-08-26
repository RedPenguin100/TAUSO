import pytest
from tests.complete.conftest import get_n_jobs

from tauso.populate.populate_fold import populate_mfe_features


@pytest.mark.parametrize("mini_sampled_data", [1000], indirect=True)
def test_mfe(mini_sampled_data, gene_to_data, dataframe_regression):
    data, feature_names = populate_mfe_features(mini_sampled_data.copy(), gene_to_data, n_jobs=get_n_jobs())
    dataframe_regression.check(data[["index_oligo"] + feature_names])


@pytest.mark.parametrize("mini_sampled_data", [1000], indirect=True)
def test_mfe_parallel_matches_baseline(mini_sampled_data, gene_to_data, dataframe_regression):
    """n_jobs=4 output must be numerically identical to the single-worker baseline."""
    data, feature_names = populate_mfe_features(mini_sampled_data.copy(), gene_to_data, n_jobs=4)
    # Reuse the existing n_jobs=1 baseline — same results expected
    dataframe_regression.check(data[["index_oligo"] + feature_names], basename="test_mfe_1000_")
