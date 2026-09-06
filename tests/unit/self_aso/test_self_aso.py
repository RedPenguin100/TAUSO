"""The oligo folded on itself, paired with a copy of itself, and overlapped end to end.

These pin behaviour that must hold whatever the weights say -- a hairpin needs a loop,
chemistry changes the answer, a mismatch does not end the whole-overlap walk. The numbers
themselves are pinned in tests/complete/regression.
"""

import numpy as np
import pandas as pd
import pytest

from tauso.features.self_aso.self_aso import (
    FEATURE_NAMES,
    calculate_self_aso,
    encode,
)
from tauso.features.self_aso.weights import (
    DINUCLEOTIDE_ENERGY,
    MISMATCH_ENERGY,
    MODIFICATION_BONUS,
    TANDEM_MISMATCH_ENERGY,
)
from tauso.populate.populate_self_aso import (
    RNA_FOLD_FEATURE,
    populate_self_aso_features,
    self_aso_feature_name,
    self_aso_feature_names,
)

MOE_GAPMER = "MMMMMddddddddddMMMMM"
CET_GAPMER = "CCCddddddddddddddCCC"


def score(sequence, pattern):
    return {k: float(v[0]) for k, v in calculate_self_aso([sequence], [pattern]).items()}


# --- what must hold whatever the weights say ---------------------------------------


def test_only_a_no_dimer():
    seq = "AAAAAAAAAAAAAAAAAAAA"
    out = score(seq, "d" * len(seq))

    assert out["hp_dg"] == 0.0
    assert out["homodimer_perfect_dg"] == 0.0
    assert out["hp_stem"] == 0.0


def test_perfect_complementarity_dimer():
    seq = "GCGCGCGCGCGCGCGCGCGC"
    out = score(seq, "d" * len(seq))

    assert out["homodimer_perfect_dg"] < -10.0
    assert out["homodimer_perfect_helix"] == 20
    # every position pairs, so the overlap walk crosses no mismatch
    assert out["homodimer_mismatched_nmm"] == 0.0
    assert out["homodimer_mismatched_len"] == 20.0


def test_hairpin_sanity():
    seq = "GGGGGGGGTTTTCCCCCCCC"

    out = score(seq, "d" * len(seq))
    assert out["hp_dg"] < 0.0
    assert out["hp_stem"] > 1


def test_hairpin_step_too_short():
    seq = "GGGGCCCC"
    assert score(seq, "d" * len(seq))["hp_dg"] == 0.0


def test_dimer_minus_hairpin():
    """Says which structure wins; a difference, not an independent search."""
    seq = "GGGGGGGGTTTTCCCCCCCC"
    out = score(seq, "d" * len(seq))

    assert out["homodimer_perfect_minus_hp"] == pytest.approx(out["homodimer_perfect_dg"] - out["hp_dg"])


def test_per_nucleotide_features():
    seq = "AAAAGCATCACTTGATCCTG"
    out = score(seq, MOE_GAPMER)

    assert out["homodimer_perfect_dg_per_nt"] == pytest.approx(out["homodimer_perfect_dg"] / len(seq))
    assert out["hp_dg_per_nt"] == pytest.approx(out["hp_dg"] / len(seq))


STEM_SEQ = "AAGGGGGCTTTTGCCCCCAA"
"""Folds to a six-pair stem: (7,12) closing, out to (2,17). Loop TTTT."""


@pytest.mark.parametrize(
    ("label", "pattern", "expected"),
    [
        ("nothing modified", "d" * 20, 0),
        ("only the closing pair", "d" * 7 + "M" + "d" * 4 + "M" + "d" * 7, 2),
        ("only the outermost pair", "ddM" + "d" * 14 + "Mdd", 2),
        ("both wings", "MMMMMddddddddddMMMMM", 6),
        ("one strand of the stem", "MMMMMMMMMMdddddddddd", 6),
    ],
)
def test_stem_modpairs(label, pattern, expected):
    """Both strands of each stem pair count, the closing and outermost pairs included.

    `one strand of the stem` is what separates counting both strands from counting one of
    them twice: the stem pairs modified nucleotides against deoxy ones throughout, so half
    of it is modified. The symmetric patterns cannot tell the two apart.
    """
    out = score(STEM_SEQ, pattern)

    assert out["hp_stem"] == 6, label
    assert out["hp_modpairs"] == expected, label


def test_stem_modpairs_ceiling():
    """A stem of n pairs spans 2n nucleotides, so that is the ceiling."""
    for pattern in ("d" * 20, MOE_GAPMER, "M" * 10 + "d" * 10, "C" * 20):
        out = score(STEM_SEQ, pattern)

        assert out["hp_modpairs"] <= 2 * out["hp_stem"], pattern


def test_stem_excludes_tails():
    """Unpaired tails do not count; a stem pair at the sequence end does.

    Both sequences carry 2'-MOE at their first and last position only. The first folds to a
    stem that stops short of them, so they are dangling ends. The second folds all the way
    out, so the same two positions are stem.
    """
    stops_short = "AAGGGGGCTTTTGCCCCCAA"
    reaches_ends = "GGGGGGGGTTTTCCCCCCCC"
    ends_only = "M" + "d" * 18 + "M"

    assert score(stops_short, ends_only)["hp_modpairs"] == 0
    assert score(reaches_ends, ends_only)["hp_modpairs"] == 2
    # fully modified, the count is exactly the stem's nucleotides, so no tail leaks in
    for seq in (stops_short, reaches_ends):
        out = score(seq, "M" * len(seq))

        assert out["hp_modpairs"] == 2 * out["hp_stem"]


# --- chemistry ---------------------------------------------------------------------


def test_cet_encoded_from_C():
    """`chemical_pattern` writes cEt as `C`. Reading it as deoxy silently mis-scores a third
    of the dataset, so the encoder is pinned here rather than trusted."""
    seq = "ACGTACGTACGTACGTACGT"
    _, _, cet = encode([seq], [CET_GAPMER])
    _, _, deoxy = encode([seq], ["d" * len(seq)])

    assert cet[0][0] == 2, "a C in the pattern must be cEt"
    assert deoxy[0][0] == 0
    assert not np.array_equal(cet[0], deoxy[0])


def test_unknown_pattern_letter_is_deoxy():
    seq = "ACGTACGTACGTACGTACGT"
    _, _, sugars = encode([seq], ["X" + "d" * (len(seq) - 1)])

    assert sugars[0][0] == 0


def test_chemistry_changes_score():
    seq = "AAAAGCATCACTTGATCCTG"
    moe = score(seq, MOE_GAPMER)
    deoxy = score(seq, "d" * len(seq))

    assert moe["homodimer_perfect_dg"] != deoxy["homodimer_perfect_dg"]
    assert moe["homodimer_perfect_modpairs"] > deoxy["homodimer_perfect_modpairs"] == 0


# --- the whole-overlap search ------------------------------------------------------


def test_mismatch_does_not_end_overlap():
    """The perfect search stops at the first mismatch; this one charges it and goes on.
    That is the whole point of the second search, so the two must disagree here."""
    seq = "AAAAGCATCACTTGATCCTG"
    out = score(seq, MOE_GAPMER)

    assert out["homodimer_mismatched_nmm"] > 0
    assert out["homodimer_mismatched_len"] > out["homodimer_perfect_helix"]


def test_full_overlap_register():
    seq = "AAAAGCATCACTTGATCCTG"
    out = score(seq, MOE_GAPMER)

    assert out["homodimer_mismatched_full_nmm"] >= out["homodimer_mismatched_nmm"], (
        "the fully overlapped register cannot cross fewer mismatches"
    )
    assert out["homodimer_mismatched_dg"] <= out["homodimer_mismatched_full_dg"], (
        "the best register cannot be worse than one of them"
    )


def test_every_mismatched_position_is_counted():
    """Nothing pairs, so every position of the fully overlapped register is a mismatch.

    The walk charges each mismatch as it steps onto it, which skips the position it starts
    on; that one is counted separately.
    """
    seq = "AAAAAAAAAAAAAAAAAAAA"
    out = score(seq, "d" * len(seq))

    assert out["homodimer_mismatched_full_nmm"] == len(seq)
    assert out["homodimer_mismatched_nmm"] <= out["homodimer_mismatched_len"]


# --- the weight tables -------------------------------------------------------------


def test_modification_bonus_scales_with_the_sugars_a_step_spans():
    """Each modified nucleotide the step spans is worth half the per-modification bonus."""
    deoxy, moe, cet = 0, 1, 2
    plain = DINUCLEOTIDE_ENERGY[deoxy, deoxy, deoxy, deoxy]
    for sugar in (moe, cet):
        one_strand = DINUCLEOTIDE_ENERGY[deoxy, deoxy, sugar, sugar] - plain
        both_strands = DINUCLEOTIDE_ENERGY[sugar, sugar, sugar, sugar] - plain
        assert np.allclose(one_strand, -MODIFICATION_BONUS[sugar])
        assert np.allclose(both_strands, -2 * MODIFICATION_BONUS[sugar])


def test_mismatch_costs_are_penalties():
    """A mismatch costs whatever the sugars are; the stacking bonus is not credited to it."""
    assert (MISMATCH_ENERGY > 0).all()
    assert (TANDEM_MISMATCH_ENERGY > 0).all()
    assert len(np.unique(MISMATCH_ENERGY)) == 1


def test_every_step_has_a_weight():
    """Whatever the fit covers, no sugar combination is left without a number."""
    assert not np.isnan(DINUCLEOTIDE_ENERGY).any()
    assert not np.isnan(MISMATCH_ENERGY).any()


# --- the populate layer ------------------------------------------------------------


def test_populate_columns():
    df = pd.DataFrame(
        {
            "aso_sequence": ["GCGCGCGCGCGCGCGCGCGC", "AAAAGCATCACTTGATCCTG"],
            "chemical_pattern": [MOE_GAPMER, CET_GAPMER],
        }
    )
    out, names = populate_self_aso_features(df, cpus=1)
    assert names == self_aso_feature_names()
    assert len(names) == len(FEATURE_NAMES) + 1
    assert RNA_FOLD_FEATURE in names
    assert not out[names].isna().any().any()


def test_populate_missing_input():
    with pytest.raises(ValueError, match="Missing columns"):
        populate_self_aso_features(pd.DataFrame({"aso_sequence": ["ACGT"]}), cpus=1)


def test_unknown_quantity_rejected():
    with pytest.raises(ValueError, match="unknown quantity"):
        self_aso_feature_name("not_a_quantity")
