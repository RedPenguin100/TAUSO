import pytest

from tauso.populate.populate_codon_usage import populate_cai, populate_enc, populate_tai

# CDS_WINDOWS must match what was used in conftest.py
CDS_WINDOWS = [20, 30, 40, 50, 60, 70]


@pytest.mark.parametrize("mini_sampled_data", [10000], indirect=True)
def test_tai_features_regression_10000(mini_sampled_data, ref_registry, dataframe_regression):
    data = mini_sampled_data.copy()
    processed_data, feature_cols = populate_tai(data, CDS_WINDOWS, ref_registry.copy())
    dataframe_regression.check(processed_data[feature_cols])


@pytest.mark.parametrize("mini_sampled_data", [10000], indirect=True)
def test_cai_features_regression_10000(mini_sampled_data, ref_registry, dataframe_regression):
    data = mini_sampled_data.copy()
    processed_data, feature_cols = populate_cai(data, CDS_WINDOWS, ref_registry.copy())
    dataframe_regression.check(processed_data[feature_cols])


@pytest.mark.parametrize("mini_sampled_data", [10000], indirect=True)
def test_enc_features_regression_10000(mini_sampled_data, ref_registry, dataframe_regression):
    data = mini_sampled_data.copy()
    processed_data, feature_cols = populate_enc(data, CDS_WINDOWS, ref_registry.copy())
    dataframe_regression.check(processed_data[feature_cols])
