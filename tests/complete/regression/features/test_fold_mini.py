import pytest
from tests.complete.conftest import get_n_jobs

from tauso.populate.populate_fold import populate_access_features, populate_mfe_features


@pytest.mark.parametrize("mini_sampled_data", [100], indirect=True)
def test_mfe_mini(mini_sampled_data, gene_to_data, dataframe_regression):
    data, feature_names = populate_mfe_features(
        mini_sampled_data.copy(),
        gene_to_data,
        n_jobs=get_n_jobs(),
    )
    dataframe_regression.check(data[["index_oligo"] + feature_names])


@pytest.mark.parametrize("mini_sampled_data", [1000], indirect=True)
def test_access_mini(mini_sampled_data, gene_to_data, dataframe_regression):
    data, feature_names = populate_access_features(
        mini_sampled_data.copy(),
        gene_to_data,
        n_jobs=get_n_jobs(),
    )
    dataframe_regression.check(data[["index_oligo"] + feature_names])
