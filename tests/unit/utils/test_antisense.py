import pytest

from tauso.util import get_antisense, get_antisense_rna


# --- get_antisense (DNA-alphabet output) ---


@pytest.mark.parametrize(
    "sense, expected",
    [
        ("A", "T"),
        ("T", "A"),
        ("G", "C"),
        ("C", "G"),
        ("ACGT", "ACGT"),
        ("AAAA", "TTTT"),
        ("GCGC", "GCGC"),
        ("ATCG", "CGAT"),
        ("AATTCCGG", "CCGGAATT"),
    ],
)
def test_get_antisense_dna(sense, expected):
    assert get_antisense(sense) == expected


def test_get_antisense_handles_u():
    # U is treated as A in the complement (RNA input → DNA output)
    assert get_antisense("U") == "A"
    assert get_antisense("AUCG") == "CGAT"


def test_get_antisense_empty():
    assert get_antisense("") == ""


def test_get_antisense_unknown_passthrough():
    # Unknown characters are passed through unchanged (str.translate skips unmapped chars)
    assert get_antisense("ACXGT") == "ACXGT"


def test_get_antisense_involution():
    # RC of RC is the original sequence
    seq = "ATCGATCG"
    assert get_antisense(get_antisense(seq)) == seq


# --- get_antisense_rna (RNA-alphabet output) ---


@pytest.mark.parametrize(
    "sense, expected",
    [
        ("A", "U"),
        ("T", "A"),
        ("G", "C"),
        ("C", "G"),
        ("U", "A"),
        ("ACGT", "ACGU"),
        ("AAAA", "UUUU"),
        ("ATCG", "CGAU"),
        ("AATTCCGG", "CCGGAAUU"),
    ],
)
def test_get_antisense_rna(sense, expected):
    assert get_antisense_rna(sense) == expected


def test_get_antisense_rna_handles_u_input():
    # RNA input with U instead of T
    assert get_antisense_rna("AUCG") == "CGAU"


def test_get_antisense_rna_empty():
    assert get_antisense_rna("") == ""


def test_get_antisense_rna_involution():
    # RC of RC is original (T becomes U in intermediate, then A, then U again)
    seq = "AUCGAUCG"
    assert get_antisense_rna(get_antisense_rna(seq)) == seq


def test_get_antisense_rna_matches_biopython():
    # Verify against Bio.Seq as reference for a handful of sequences
    from Bio.Seq import Seq

    seqs = ["ACGT", "ATCGATCG", "GCGCGCGC", "AAAAGGGG", "TTTTCCCC"]
    for s in seqs:
        expected = str(Seq(s).reverse_complement_rna())
        assert get_antisense_rna(s) == expected, f"mismatch for {s!r}"
