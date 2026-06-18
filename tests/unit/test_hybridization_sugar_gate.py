"""The cEt/LNA hybridization features apply only to a single-high-affinity-sugar
oligo (that sugar + DNA). Mixmers, mixed chemistries, and unknown codes -> NaN."""

import math

import pytest

from tauso.features.hybridization.hybridization_features import (
    get_cet_dna_rna_dg,
    get_cet_wing_dg,
    get_lna_dna_rna_dg,
    is_single_sugar_mod,
)

SEQ = "ACGTACGTACGTACGTACGT"


def _is_nan(x):
    return isinstance(x, float) and math.isnan(x)


def test_single_sugar_mod_predicate():
    assert is_single_sugar_mod("CCCddddddddddCCC", "C")
    assert not is_single_sugar_mod("MMCddddddddddCCM", "C")  # mixmer (also MOE)
    assert not is_single_sugar_mod("CCLddddddddddLCC", "C")  # also LNA
    assert not is_single_sugar_mod("CCXddddddddddXCC", "C")  # unknown code
    assert not is_single_sugar_mod("dddddddddddddddd", "C")  # no cEt
    assert not is_single_sugar_mod(None, "C")
    assert is_single_sugar_mod("LLLddddddddddLLL", "L")


@pytest.mark.parametrize(
    "pattern",
    [
        "MMCddddddddddCCM",  # cEt + MOE (mixmer)
        "CCLddddddddddLCC",  # cEt + LNA
        "CCXddddddddddXCC",  # cEt + unknown
        "MMMddddddddddMMM",  # no cEt
        "dddddddddddddddd",  # PS-DNA
    ],
)
def test_cet_is_nan_unless_pure_cet(pattern):
    a = SEQ[: len(pattern)]
    assert _is_nan(get_cet_dna_rna_dg(a, pattern))
    assert _is_nan(get_cet_wing_dg(a, pattern, "wing5"))


def test_cet_populated_for_pure_cet():
    p = "CCCddddddddddCCC"
    a = SEQ[: len(p)]
    assert not _is_nan(get_cet_dna_rna_dg(a, p))
    assert not _is_nan(get_cet_wing_dg(a, p, "wing5"))


def test_lna_gate_is_sole_sugar():
    a = SEQ[:16]
    assert _is_nan(get_lna_dna_rna_dg(a, "CCLddddddddddLCC"))  # LNA + cEt -> not pure LNA
    assert not _is_nan(get_lna_dna_rna_dg(a, "LLLddddddddddLLL"))  # pure LNA
