import pandas as pd
import pytest

from tauso.populate.populate_rbp import (
    populate_complexity_features,
    populate_functional_features,
    populate_rbp_affinity_features,
    populate_rbp_interaction_features,
)
from tests.complete.conftest import get_n_jobs


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
            data, rbp_map, pwm_db, gene_to_data, window_col, n_jobs=get_n_jobs()
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


@pytest.mark.parametrize("sampled_base_data", [1000], indirect=True)
def test_rbp_interaction_features_regression(
    sampled_base_data, gene_to_data, rbp_assets, expression_matrix, dataframe_regression
):
    """
    Tests the expression-weighted interaction path -- the only consumer of
    build_rbp_expression_matrix (pulled in via the expression_matrix fixture) --
    plus the complexity (expression) and functional (stabilizer/destabilizer)
    aggregations, mirroring calculator.calculate_rbp blocks A and C.
    """
    from tauso.features.rbp.rbp_annotations import rbp_role_map_strict

    data = sampled_base_data.copy()
    rbp_map, pwm_db = rbp_assets

    flank_size = 50
    window_col = f"flank_sequence_{flank_size}"

    data, interaction_features = populate_rbp_interaction_features(
        data, rbp_map, pwm_db, expression_matrix, gene_to_data, window_col, n_jobs=get_n_jobs()
    )

    data, global_features = populate_complexity_features(
        data, interaction_features, suffix=str(flank_size), type="expression"
    )

    data, functional_features = populate_functional_features(data, rbp_role_map_strict, flank_size=flank_size)

    dataframe_regression.check(
        data[interaction_features + global_features + functional_features],
        default_tolerance=dict(atol=1e-5, rtol=1e-5),
    )
