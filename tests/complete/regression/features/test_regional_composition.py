import pytest

from tauso.populate.populate_regional_composition import populate_regional_composition_features


@pytest.mark.parametrize("mini_sampled_data", [2000], indirect=True)
def test_regional_composition_features(mini_sampled_data, dataframe_regression):
    data = mini_sampled_data.copy()
    processed_data, feature_cols = populate_regional_composition_features(data)
    dataframe_regression.check(processed_data[feature_cols])
