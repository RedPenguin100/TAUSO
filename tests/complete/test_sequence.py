import pytest

from tauso.populate.populate_sequence import populate_sequence_features


@pytest.fixture
def mini_sampled_data(request, final_data):
    """Samples the fully processed DataFrame right before the test runs."""
    n_samples = getattr(request, "param", 1000)
    actual_samples = min(n_samples, len(final_data))
    return final_data.sample(n=actual_samples, random_state=42).copy()


@pytest.mark.parametrize("mini_sampled_data", [10000], indirect=True)
def test_sequence_features(mini_sampled_data, dataframe_regression):
    data = mini_sampled_data.copy()
    processed_data, feature_cols = populate_sequence_features(data)
    dataframe_regression.check(processed_data[feature_cols])
