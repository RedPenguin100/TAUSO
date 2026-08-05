import numpy as np
import pytest

from tauso.util import _norm_rna_to_dna, _to_str_seq, normalize_dna


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("ACGT", "ACGT"),
        ("acgt", "ACGT"),
        ("ACGU", "ACGT"),  # RNA U -> DNA T
        ("acgu", "ACGT"),  # lowercase u -> T (uppercase then substitute)
        (" A C\tG\nT ", "ACGT"),  # whitespace stripped
        ("", ""),
    ],
)
def test_norm_rna_to_dna(raw, expected):
    assert _norm_rna_to_dna(raw) == expected


def test_to_str_seq_matches_norm_on_strings():
    for raw in ["ACGT", "acgu", " a c g u ", "UUUU", "GgCcTtAa"]:
        assert _to_str_seq(raw) == _norm_rna_to_dna(raw)


@pytest.mark.parametrize(
    "seq_like, expected",
    [
        (["A", "C", "G", "U"], "ACGT"),
        (np.array(["a", "c", "g", "t"]), "ACGT"),
        (("U", "U", "U"), "TTT"),
    ],
)
def test_to_str_seq_coerces_sequence_like(seq_like, expected):
    assert _to_str_seq(seq_like) == expected


@pytest.mark.parametrize("raw, expected", [("acgu", "ACGT"), ("ACGT", "ACGT"), (np.str_("acgt"), "ACGT")])
def test_normalize_dna_accepts_str_like(raw, expected):
    assert normalize_dna(raw) == expected


@pytest.mark.parametrize("raw, match", [("ACG ACG", "whitespace"), ("", "empty"), ("ACGZT", "non-ACGU")])
def test_normalize_dna_rejects_bad_sequences(raw, match):
    with pytest.raises(ValueError, match=match):
        normalize_dna(raw)


@pytest.mark.parametrize("raw", [None, float("nan"), 42, True])
def test_normalize_dna_rejects_non_string_input(raw):
    """Non-strings stringify into letters ('nan' -> 'NAN'), so they must fail as a type error
    rather than be reported as a bad base."""
    with pytest.raises(TypeError, match="must be a string"):
        normalize_dna(raw)
