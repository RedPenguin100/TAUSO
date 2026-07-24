"""The guard that stops candidates with failed feature computations from being scored.

Part of the model's feature set is legitimately missing for part of the training data, so the
booster learned a missing-value branch for it; the rest is always present, and a NaN there can only
be a failed computation. These tests cover the split and the detection, without needing the booster.
"""

import numpy as np
import pytest

from tauso.inference import DEFAULT_VERSION, load_finite_features
from tauso.inference.predict import MODEL_DIR, _classify_nonfinite, _find_nonfinite


def model_features():
    path = MODEL_DIR / f"tauso_score_{DEFAULT_VERSION}.features.txt"
    return [f for f in path.read_text().splitlines() if f.strip()]


def test_finite_features_are_a_proper_subset_of_the_model_features():
    features = model_features()
    finite = load_finite_features()
    assert finite  # the guard is armed
    assert finite < set(features)  # every guarded name is a model feature, and some are exempt


def test_conditionally_defined_features_are_exempt():
    """Codon scores (CDS windows only) and cEt hybridization terms (cEt chemistries only) are NaN
    for most of the training data, so the booster handles them and the guard must not fire on them."""
    finite = load_finite_features()
    assert "cai_score_20" not in finite
    assert "hybr_cet_dna_rna_dg" not in finite


@pytest.fixture
def split():
    """(features, one guarded name, one exempt name) for building matrices to check."""
    features = model_features()
    finite = load_finite_features()
    guarded = next(f for f in features if f in finite)
    exempt = next(f for f in features if f not in finite)
    return features, guarded, exempt


def test_finite_matrix_is_clean(split):
    features, _, _ = split
    X = np.zeros((3, len(features)))
    assert _find_nonfinite(X, features, DEFAULT_VERSION) == []


def test_nan_in_a_guarded_feature_is_reported(split):
    features, guarded, _ = split
    X = np.zeros((3, len(features)))
    X[[0, 2], features.index(guarded)] = np.nan
    assert _find_nonfinite(X, features, DEFAULT_VERSION) == [(guarded, 2)]


def test_nan_in_an_exempt_feature_is_ignored(split):
    features, _, exempt = split
    X = np.zeros((3, len(features)))
    X[:, features.index(exempt)] = np.nan
    assert _find_nonfinite(X, features, DEFAULT_VERSION) == []


def test_infinity_is_reported_even_for_an_exempt_feature(split):
    """No feature is ever infinite in training, so an infinity is out of distribution regardless of
    whether the feature is allowed to be missing."""
    features, _, exempt = split
    X = np.zeros((3, len(features)))
    X[0, features.index(exempt)] = np.inf
    X[1, features.index(exempt)] = -np.inf
    assert _find_nonfinite(X, features, DEFAULT_VERSION) == [(exempt, 2)]


def test_offenders_are_ordered_worst_first(split):
    features, guarded, exempt = split
    other = next(f for f in features if f in load_finite_features() and f != guarded)
    X = np.zeros((4, len(features)))
    X[[0], features.index(guarded)] = np.nan
    X[[0, 1, 2], features.index(other)] = np.nan
    X[[0, 1], features.index(exempt)] = np.inf
    assert _find_nonfinite(X, features, DEFAULT_VERSION) == [(other, 3), (exempt, 2), (guarded, 1)]


def test_partial_batch_is_a_failure_and_whole_batch_is_a_gap():
    """A feature that computed for some candidates but not others failed for those; one missing for
    every candidate is unavailable for the target."""
    offenders = [("gene_level_annotation", 10), ("failed_fold", 3)]
    failures, gaps = _classify_nonfinite(offenders, n_candidates=10)
    assert failures == [("failed_fold", 3)]
    assert gaps == [("gene_level_annotation", 10)]


def test_single_candidate_reads_as_a_gap():
    """With one candidate there is nothing to be inconsistent with, so nothing is a failure."""
    failures, gaps = _classify_nonfinite([("anything", 1)], n_candidates=1)
    assert failures == []
    assert gaps == [("anything", 1)]


def test_clean_batch_classifies_to_nothing():
    assert _classify_nonfinite([], n_candidates=10) == ([], [])
