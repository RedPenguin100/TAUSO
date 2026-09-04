import pytest

from tauso.populate.populate_self_aso import populate_self_aso_features


@pytest.mark.parametrize("mini_sampled_data", [10000], indirect=True)
def test_self_aso_features(mini_sampled_data, dataframe_regression):
    data = mini_sampled_data.copy()
    processed_data, feature_cols = populate_self_aso_features(data)
    dataframe_regression.check(processed_data[feature_cols])
