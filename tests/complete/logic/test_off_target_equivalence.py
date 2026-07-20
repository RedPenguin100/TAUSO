import numpy as np
import pandas as pd
import pytest
from tests.complete.conftest import get_n_jobs

from tauso.features.hybridization.off_target.off_target_specific_gene import (
    add_on_target_site_features,
    on_target_total_hybridization,
)
from tauso.features.hybridization.off_target.on_target_multiplicity import (
    on_target_log_number_of_sites,
)

CUTOFFS = [800, 1000, 1200]


@pytest.mark.parametrize("mini_structure_data", [1000], indirect=True)
def test_add_on_target_site_features_matches_separate(mini_structure_data, gene_to_data):
    """The single-scan add_on_target_site_features emits the same total-hybridization and
    log-number-of-sites columns as running the two per-feature entry points separately."""
    n_jobs = get_n_jobs()
    both, both_names = add_on_target_site_features(
        mini_structure_data.copy(), gene_to_data, cutoffs=CUTOFFS, n_jobs=n_jobs
    )
    tot, tot_names = on_target_total_hybridization(
        mini_structure_data.copy(), gene_to_data, cutoffs=CUTOFFS, n_jobs=n_jobs
    )
    log, log_names = on_target_log_number_of_sites(
        mini_structure_data.copy(), gene_to_data, cutoffs=CUTOFFS, n_jobs=n_jobs
    )
    for name in tot_names:
        assert name in both_names
        np.testing.assert_allclose(both[name].to_numpy(), tot[name].to_numpy(), rtol=1e-4, atol=1e-4)
    for name in log_names:
        assert name in both_names
        np.testing.assert_allclose(both[name].to_numpy(), log[name].to_numpy(), rtol=1e-4, atol=1e-4)


@pytest.mark.parametrize("mini_structure_data", [300], indirect=True)
def test_on_target_multi_cutoff_matches_per_cutoff(mini_structure_data, gene_to_data):
    """on_target deriving all cutoffs from one loose pass equals running each separately
    (within FP rounding of the order-dependent exp-sum)."""
    cutoffs = [800, 1000, 1200]
    multi, _ = on_target_total_hybridization(
        mini_structure_data.copy(), gene_to_data, cutoffs=cutoffs, n_jobs=get_n_jobs()
    )
    for c in cutoffs:
        single, sfeats = on_target_total_hybridization(
            mini_structure_data.copy(), gene_to_data, cutoffs=[c], n_jobs=get_n_jobs()
        )
        pd.testing.assert_series_equal(multi[sfeats[0]], single[sfeats[0]])
