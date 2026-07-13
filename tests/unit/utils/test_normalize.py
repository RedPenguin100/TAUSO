import numpy as np
import pytest

from tauso.util import _norm_rna_to_dna, _to_str_seq


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
