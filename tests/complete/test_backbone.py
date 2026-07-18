import pytest

from tauso.populate.populate_backbone import populate_backbone_features


@pytest.mark.parametrize("sampled_base_data", [10000], indirect=True)
def test_backbone_features(sampled_base_data, dataframe_regression):
    data = sampled_base_data.copy()
    processed_data, feature_cols = populate_backbone_features(data)
    dataframe_regression.check(processed_data[feature_cols])
