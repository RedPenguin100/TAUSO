import pytest

from tauso.features.context.mrna_halflife import (
    populate_mrna_halflife_features,
)


@pytest.mark.parametrize("mini_sampled_data", [10000], indirect=True)
def test_half_life(mini_sampled_data, half_life_provider, dataframe_regression):
    data = mini_sampled_data.copy()

    processed_data, feature_cols = populate_mrna_halflife_features(data, half_life_provider)

    dataframe_regression.check(processed_data[feature_cols])
