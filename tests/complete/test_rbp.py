import pandas as pd
import pytest

from tauso.populate.populate_rbp import (
    populate_complexity_features,
    populate_rbp_affinity_features,
)


@pytest.fixture
def sampled_base_data(request, final_data):
    """Samples the base data right before the test runs."""
    n_samples = getattr(request, "param", 1000)
    actual_samples = min(n_samples, len(final_data))
    return final_data.sample(n=actual_samples, random_state=42).copy()


@pytest.mark.parametrize("sampled_base_data", [1000], indirect=True)
def test_rbp_complexity_features_regression(
    sampled_base_data, gene_to_data, rbp_assets, dataframe_regression, data_regression
):
    """
    Tests the RBP affinity and complexity feature population.
    Uses dataframe_regression for the main df, and data_regression for the feature dicts/dfs.
    """
    data = sampled_base_data.copy()
    rbp_map, pwm_db = rbp_assets

    for flank_size in [50]:
        window_col = f"flank_sequence_{flank_size}"

        data, individual_features = populate_rbp_affinity_features(
            data, rbp_map, pwm_db, gene_to_data, window_col, n_jobs=2
        )

        data, global_features = populate_complexity_features(
            data, individual_features, suffix=str(flank_size), type="generic"
        )

    # DataFrame Regression
    dataframe_regression.check(
        data[individual_features + global_features],
        default_tolerance=dict(atol=1e-5, rtol=1e-5),
    )

    # Data Regression for complex payload
    feature_payload = {
        "individual_features": individual_features.to_dict("records")
        if isinstance(individual_features, pd.DataFrame)
        else individual_features,
        "global_features": global_features.to_dict("records")
        if isinstance(global_features, pd.DataFrame)
        else global_features,
    }

    data_regression.check(feature_payload)
