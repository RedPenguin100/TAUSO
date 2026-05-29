"""Mini fold/access regression test at 100 samples.

Wired to the module-level defaults of populate_fold so that re-parametrizing
the defaults forces a baseline regen. Also reports wall-clock so we can
compare runtime across parameter sets.
"""
import time
from collections import defaultdict

import pytest

from tauso.populate.populate_fold import (
    DEFAULT_SENSE_CONFIGURATION,
    populate_mfe_features,
    populate_sense_accessibility_multi,
)
from tests.complete.conftest import get_n_jobs


@pytest.mark.parametrize("mini_sampled_data", [100], indirect=True)
def test_mfe_mini(mini_sampled_data, gene_to_data, dataframe_regression):
    t0 = time.perf_counter()
    data, feature_names = populate_mfe_features(
        mini_sampled_data.copy(),
        gene_to_data,
        n_jobs=get_n_jobs(),
    )
    elapsed = time.perf_counter() - t0
    print(
        f"\n[MFE mini] n_jobs={get_n_jobs()} "
        f"n_rows={len(mini_sampled_data)} "
        f"features={len(feature_names)} "
        f"elapsed={elapsed:.2f}s "
        f"per_row={elapsed / len(mini_sampled_data) * 1000:.1f}ms"
    )
    dataframe_regression.check(data[["index_oligo"] + feature_names])


@pytest.mark.parametrize("mini_sampled_data", [100], indirect=True)
def test_access_mini(mini_sampled_data, gene_to_data, dataframe_regression):
    # populate_sense_accessibility_multi requires single-flank configs per call;
    # group DEFAULT_SENSE_CONFIGURATION by flank and check baselines per global
    # config index.
    groups = defaultdict(list)
    for idx, cfg in enumerate(DEFAULT_SENSE_CONFIGURATION):
        groups[cfg["flank"]].append((idx, cfg))

    t0 = time.perf_counter()
    all_feature_names = []
    data = mini_sampled_data
    for flank, indexed_cfgs in groups.items():
        idxs, configs = zip(*indexed_cfgs)
        data, feature_names = populate_sense_accessibility_multi(
            data,
            gene_to_data,
            configs=list(configs),
            batch_size=100,
            n_jobs=get_n_jobs(),
        )
        for idx, feature_name in zip(idxs, feature_names):
            dataframe_regression.check(
                data[["index_oligo", feature_name]],
                basename=f"test_access_mini_config{idx}",
            )
            all_feature_names.append(feature_name)
    elapsed = time.perf_counter() - t0
    print(
        f"\n[Access mini] n_jobs={get_n_jobs()} "
        f"n_rows={len(mini_sampled_data)} "
        f"configs={len(DEFAULT_SENSE_CONFIGURATION)} "
        f"features={len(all_feature_names)} "
        f"elapsed={elapsed:.2f}s "
        f"per_row={elapsed / len(mini_sampled_data) * 1000:.1f}ms"
    )
