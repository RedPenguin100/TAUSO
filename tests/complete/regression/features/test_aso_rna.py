import pytest

from tauso.populate.populate_aso_rna import populate_aso_rna_features


@pytest.mark.parametrize("mini_sampled_data", [2000], indirect=True)
def test_aso_rna_features(mini_sampled_data, dataframe_regression):
    data = mini_sampled_data.copy()
    processed_data, feature_cols = populate_aso_rna_features(data)
    dataframe_regression.check(processed_data[feature_cols])
