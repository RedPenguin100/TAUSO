import numpy as np
import pandas as pd
import pytest

from tauso.inference import DEFAULT_VERSION, load_finite_features, load_model, predict


@pytest.fixture(scope="module")
def features():
    _, names = load_model(DEFAULT_VERSION)
    return names


@pytest.fixture(scope="module")
def guarded_and_unguarded(features):
    guarded = load_finite_features(DEFAULT_VERSION)
    return next(f for f in features if f in guarded), next(f for f in features if f not in guarded)


@pytest.fixture
def frame(features):
    return pd.DataFrame(np.zeros((4, len(features))), columns=features)


def test_clean_frame_scores(frame):
    assert len(predict(frame)) == 4


def test_partial_failure_raises(frame, guarded_and_unguarded):
    guarded, _ = guarded_and_unguarded
    frame.loc[0, guarded] = np.nan
    with pytest.raises(ValueError, match="failed to compute for part of"):
        predict(frame)


def test_partial_failure_is_a_warning_when_not_strict(frame, guarded_and_unguarded, caplog):
    guarded, _ = guarded_and_unguarded
    frame.loc[0, guarded] = np.nan
    scores = predict(frame, strict=False)
    assert len(scores) == 4
    assert "failed to compute for part of" in caplog.text


def test_whole_batch_gap_warns_but_scores(frame, guarded_and_unguarded, caplog):
    """A feature missing for every candidate is unavailable for the target, not broken."""
    guarded, _ = guarded_and_unguarded
    frame[guarded] = np.nan
    scores = predict(frame)
    assert len(scores) == 4
    assert "non-finite for all 4 candidates" in caplog.text


def test_nan_in_an_unguarded_feature_is_silent(frame, guarded_and_unguarded, caplog):
    _, unguarded = guarded_and_unguarded
    frame.loc[0, unguarded] = np.nan
    assert len(predict(frame)) == 4
    assert caplog.text == ""


def test_infinity_raises_even_when_unguarded(frame, guarded_and_unguarded):
    _, unguarded = guarded_and_unguarded
    frame.loc[0, unguarded] = np.inf
    with pytest.raises(ValueError, match="failed to compute for part of"):
        predict(frame)


def test_missing_column_still_raises_first(frame, features):
    with pytest.raises(ValueError, match="missing"):
        predict(frame.drop(columns=[features[0]]))
