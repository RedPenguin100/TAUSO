import logging

import pytest

from tauso.populate.populate_fold import (
    populate_mfe_features,
    populate_sense_accessibility_batch,
)
from tests.complete.conftest import get_n_jobs


ACCESS_CONFIGURATIONS = [
    {"flank": 120, "access": 20, "seeds": [13]},
    {"flank": 120, "access": 13, "seeds": [4, 6, 8]},
    {"flank": 120, "access": 20, "seeds": [4, 6, 8]},
    {"flank": 120, "access": 13, "seeds": [13, 26, 39]},
]


@pytest.mark.parametrize("mini_sampled_data", [1000], indirect=True)
def test_access(mini_sampled_data, gene_to_data, dataframe_regression):
    for idx, config in enumerate(ACCESS_CONFIGURATIONS):
        data, feature_name = populate_sense_accessibility_batch(
            mini_sampled_data,
            gene_to_data,
            batch_size=250,
            flank_size=config["flank"],
            access_size=config["access"],
            seed_sizes=config["seeds"],
            n_jobs=get_n_jobs(),
        )
        dataframe_regression.check(
            data[["index_oligo", feature_name]],
            basename=f"test_access_config{idx}_1000_",
        )


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

