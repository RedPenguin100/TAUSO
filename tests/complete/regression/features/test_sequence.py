import pytest

from tauso.populate.populate_sequence import populate_sequence_features


@pytest.mark.parametrize("mini_sampled_data", [10000], indirect=True)
def test_sequence_features(mini_sampled_data, dataframe_regression):
    data = mini_sampled_data.copy()
    processed_data, feature_cols = populate_sequence_features(data)
    dataframe_regression.check(processed_data[feature_cols])
