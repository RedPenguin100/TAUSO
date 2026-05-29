import logging
from collections import defaultdict

import pytest

from tauso.populate.populate_fold import (
    DEFAULT_SENSE_CONFIGURATION,
    populate_mfe_features,
    populate_sense_accessibility_multi,
)
from tests.complete.conftest import get_n_jobs


@pytest.mark.parametrize("mini_sampled_data", [1000], indirect=True)
def test_access(mini_sampled_data, gene_to_data, dataframe_regression):
    # populate_sense_accessibility_multi requires all configs in one call to share
    # `flank` (one raccess pass per batch). DEFAULT_SENSE_CONFIGURATION spans
    # several flanks, so group + run once per flank; baselines stay one-per-config
    # via the global config index used as the basename.
    groups = defaultdict(list)
    for idx, cfg in enumerate(DEFAULT_SENSE_CONFIGURATION):
        groups[cfg["flank"]].append((idx, cfg))

    data = mini_sampled_data
    for flank, indexed_cfgs in groups.items():
        idxs, configs = zip(*indexed_cfgs)
        data, feature_names = populate_sense_accessibility_multi(
            data,
            gene_to_data,
            configs=list(configs),
            batch_size=250,
            n_jobs=get_n_jobs(),
        )
        for idx, feature_name in zip(idxs, feature_names):
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

