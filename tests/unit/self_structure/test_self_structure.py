"""The ASO folded on itself, and paired with a second copy of itself.

The weights are a fitted table, so these pin structure rather than values: that a hairpin
needs a loop, that chemistry changes the answer, and that the reported energies are the best
structure found rather than an arbitrary one.
"""

import numpy as np
import pytest

from tauso.features.self_structure.self_structure import (
    FEATURE_NAMES,
    calculate_self_structure,
    encode,
)
from tauso.populate.populate_self_structure import self_structure_feature_name

DEOXY_20 = "D" * 20


def score(sequence, pattern):
    return {k: v[0] for k, v in calculate_self_structure([sequence], [pattern]).items()}


def test_a_sequence_that_cannot_pair_has_no_structure():
    out = score("AAAAAAAAAAAAAAAAAAAA", DEOXY_20)
    assert out["hp_dg"] == 0.0
    assert out["dim_dg"] == 0.0
    assert out["hp_stem"] == 0.0


def test_a_self_complementary_oligo_dimerises_strongly():
    out = score("GCGCGCGCGCGCGCGCGCGC", DEOXY_20)
    assert out["dim_dg"] < -10.0
    assert out["dim_helix"] > 10


def test_a_stem_loop_folds():
    stem_loop = score("GGGGGGGGTTTTCCCCCCCC", DEOXY_20)
    assert stem_loop["hp_dg"] < 0.0
    assert stem_loop["hp_stem"] > 1


def test_a_short_stem_does_not_pay_for_its_loop():
    """A hairpin closes a loop; a stem too short to repay that penalty is not reported."""
    assert score("GGGGCCCC", "D" * 8)["hp_dg"] == 0.0
    # The same arms with a longer stem clear it.
    assert score("GGGGGGGGTTTTCCCCCCCC", DEOXY_20)["hp_dg"] < 0.0
    # A dimer pays no loop penalty, so the short one still pairs.
    assert score("GGGGCCCC", "D" * 8)["dim_dg"] < 0.0


def test_energies_are_the_best_structure_so_never_positive():
    seqs = ["GCATTGGTATTCAGTGTGAT", "CCTTCCCTGAAGGTTCCTCC", "GGGGGGGGTTTTCCCCCCCC"]
    out = calculate_self_structure(seqs, [DEOXY_20] * 3)
    assert (out["hp_dg"] <= 0).all()
    assert (out["dim_dg"] <= 0).all()


def test_chemistry_changes_the_energy():
    """2'-MOE and cEt carry their own stacking rather than being read as deoxy."""
    seq = "GGGGGGGGTTTTCCCCCCCC"
    deoxy = score(seq, DEOXY_20)["hp_dg"]
    moe = score(seq, "MMMMM" + "D" * 10 + "MMMMM")["hp_dg"]
    cet = score(seq, "CCCCC" + "D" * 10 + "CCCCC")["hp_dg"]
    assert deoxy != moe
    assert moe != cet


def test_dim_minus_hp_is_the_difference():
    out = score("GGGGGGGGTTTTCCCCCCCC", DEOXY_20)
    assert out["dim_minus_hp"] == pytest.approx(out["dim_dg"] - out["hp_dg"])


def test_per_nucleotide_divides_by_length():
    seq = "GGGGGGGGTTTTCCCCCCCC"
    out = score(seq, DEOXY_20)
    assert out["hp_dg_per_nt"] == pytest.approx(out["hp_dg"] / len(seq))
    assert out["dim_dg_per_nt"] == pytest.approx(out["dim_dg"] / len(seq))


def test_a_base_outside_acgt_is_nan_not_zero():
    """Zero would read as a free oligo; an unstackable base is unknown, not free."""
    out = calculate_self_structure(["ACGNACGT"], ["D" * 8])
    for name in FEATURE_NAMES:
        assert np.isnan(out[name]).all()


def test_length_is_not_capped():
    """The buffer follows the data, so a longer oligo scores rather than silently dropping."""
    long_stem_loop = "GGGGGGGGGGGGTTTTCCCCCCCCCCCC"
    assert len(long_stem_loop) > 21
    assert score(long_stem_loop, "D" * len(long_stem_loop))["hp_dg"] < 0.0


def test_encode_marks_moe_and_cet_and_leaves_the_rest_deoxy():
    _, lengths, sugars = encode(["ACGT"], ["MCDx"])
    assert lengths[0] == 4
    assert list(sugars[0][:4]) == [1, 2, 0, 0]


def test_feature_names_carry_the_weight_set():
    assert self_structure_feature_name("hp_dg") == "md3_hp_dg"
    assert self_structure_feature_name("hp_dg", "md4") == "md4_hp_dg"
    with pytest.raises(ValueError, match="unknown quantity"):
        self_structure_feature_name("hp_enthalpy")
