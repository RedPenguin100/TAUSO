import pytest

from tauso.populate.populate_fold import (
    populate_mfe_features,
    populate_sense_accessibility_batch,
)

ACCESS_CONFIGURATIONS = [
    {"flank": 120, "access": 20, "seeds": [13]},
    {"flank": 120, "access": 13, "seeds": [4, 6, 8]},
    {"flank": 120, "access": 20, "seeds": [4, 6, 8]},
    {"flank": 120, "access": 13, "seeds": [13, 26, 39]},
]


@pytest.mark.parametrize("mini_sampled_data", [1000], indirect=True)
@pytest.mark.parametrize("config", ACCESS_CONFIGURATIONS)
def test_access(mini_sampled_data, gene_to_data, config, dataframe_regression):
    data, feature_name = populate_sense_accessibility_batch(
        mini_sampled_data,
        gene_to_data,
        batch_size=1000,
        flank_size=config["flank"],
        access_size=config["access"],
        seed_sizes=config["seeds"],
        n_jobs=1,
    )
    dataframe_regression.check(data[["index_oligo", feature_name]])


@pytest.mark.parametrize("mini_sampled_data", [1000], indirect=True)
def test_mfe(mini_sampled_data, gene_to_data, dataframe_regression):
    # populate_mfe_features mutates its input — copy to protect the session fixture
    data, feature_names = populate_mfe_features(mini_sampled_data.copy(), gene_to_data)
    dataframe_regression.check(data[["index_oligo"] + feature_names])
