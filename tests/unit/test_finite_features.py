import pytest

from tauso.inference import DEFAULT_VERSION, load_finite_features, load_model


def test_finite_features_are_model_features():
    finite = load_finite_features()
    _, model_features = load_model(DEFAULT_VERSION)
    assert finite, "the finite-feature list is empty"
    assert finite <= set(model_features)


def test_some_model_features_are_left_unguarded():
    """The conditionally-defined features (codon scores, cEt terms, junction distances) are missing
    in training by design, so guarding every one would reject legitimate input."""
    _, model_features = load_model(DEFAULT_VERSION)
    unguarded = set(model_features) - load_finite_features()
    assert unguarded
    assert any(f.startswith(("cai_score", "tai_score", "enc_score", "hybr_cet")) for f in unguarded)


def test_unknown_version_raises():
    with pytest.raises(FileNotFoundError):
        load_finite_features("v_does_not_exist")
