"""What the bundled model asks for that the feature pipeline no longer produces.

`predict` selects the model's feature list out of the frame and raises on anything missing, so
a retired feature makes the bundled model unscorable rather than slightly wrong. Two rounds of
feature work have retired columns v1 was trained on, and it has not been retrained since.

This records exactly which ones, so the gap cannot grow unnoticed and whoever retrains knows
what the new model has to do without.
"""

import pytest

from tauso.inference import DEFAULT_VERSION, load_finite_features, load_model
from tauso.populate.populate_self_aso import self_aso_feature_names
from tauso.populate.populate_sequence import FEATURE_SPECS

RETIRED_SELF_STRUCTURE = frozenset(
    {
        "seq_self_energy",
        "seq_hairpin_score",
        "seq_hairpin_dG_energy",
        "seq_hairpin_tm",
        "interaction_internal_fold_rna_gymnosis",
    }
)
"""Scored the oligo's own structure before the self-ASO family replaced them."""

RETIRED_ACCESS_PREFIX = "fold_access_"
"""The accessibility family that went with raccess, replaced by the ViennaRNA one."""


def _retired(model_features):
    return {f for f in model_features if f in RETIRED_SELF_STRUCTURE or f.startswith(RETIRED_ACCESS_PREFIX)}


def test_the_bundled_model_is_waiting_on_a_retrain():
    """v1 cannot be scored until then. Change this number only by retraining or retiring more."""
    _, model_features = load_model(DEFAULT_VERSION)
    retired = _retired(model_features)

    assert len(retired) == 15, sorted(retired)
    assert sum(f.startswith(RETIRED_ACCESS_PREFIX) for f in retired) == 10
    assert RETIRED_SELF_STRUCTURE <= retired


def test_the_finite_guard_still_names_retired_features():
    """`finite.txt` rejects a NaN in any feature it lists, so it goes stale the same way."""
    finite = load_finite_features()

    assert RETIRED_SELF_STRUCTURE <= finite


@pytest.mark.parametrize("name", sorted(RETIRED_SELF_STRUCTURE))
def test_a_retired_feature_does_not_come_back(name):
    """Guards the retirement itself: nothing should start producing these again by accident."""
    produced = {n for n, _ in FEATURE_SPECS} | set(self_aso_feature_names())

    assert name not in produced
