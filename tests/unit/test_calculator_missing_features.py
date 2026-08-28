import pandas as pd
import pytest

from tauso.data.consts import CANONICAL_GENE_NAME
from tauso.populate.calculators.calculator import Calculator


@pytest.fixture
def calculator_factory():
    def make(**kwargs):
        data = pd.DataFrame({CANONICAL_GENE_NAME: ["MALAT1", "ACTB"]})
        return Calculator(data=data, data_version="test", **kwargs)

    return make


def test_everything_is_missing_without_a_feature_dir(calculator_factory):
    calc = calculator_factory(overwrite=False, get_feature_dir=None)

    assert calc._get_missing_features(["feat_a", "feat_b"]) == ["feat_a", "feat_b"]


def test_saved_features_are_not_missing(calculator_factory, tmp_path):
    calc = calculator_factory(overwrite=False, get_feature_dir=lambda version: str(tmp_path))
    calc.data["feat_a"] = [1.0, 2.0]
    calc.data["index_test"] = [0, 1]
    calc._save_calculated_feature(feature_name="feat_a")

    assert calc._get_missing_features(["feat_a", "feat_b"]) == ["feat_b"]


def test_overwrite_reports_everything_missing(calculator_factory, tmp_path):
    calc = calculator_factory(overwrite=True, get_feature_dir=lambda version: str(tmp_path))

    assert calc._get_missing_features(["feat_a"]) == ["feat_a"]
