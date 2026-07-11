"""Target-duplication matching: exact vs one-mismatch counting and clustering."""

import numpy as np
import pytest

from tauso.features.duplication_features import compute_duplications, distinct_matches

TARGET = "ACGUACGUACGUACGU"  # 16 nt


def _wrap(core):
    return "UU" + core + "UU"


def test_exact_single_copy():
    mrna = _wrap(TARGET)
    assert distinct_matches(mrna, TARGET, 0) == 1
    assert distinct_matches(mrna, TARGET, 1) == 1


def test_absent_target():
    mrna = "U" * 60
    assert distinct_matches(mrna, TARGET, 0) == 0
    assert distinct_matches(mrna, TARGET, 1) == 0


def test_one_mismatch_left_half():
    # substitution in the left half: only the right seed anchors the hit
    mm = "ACGUACGCACGUACGU"
    mrna = _wrap(mm)
    assert distinct_matches(mrna, TARGET, 0) == 0
    assert distinct_matches(mrna, TARGET, 1) == 1


def test_one_mismatch_right_half():
    # substitution in the right half: only the left seed anchors the hit
    mm = "ACGUACGUACGUACGC"
    mrna = _wrap(mm)
    assert distinct_matches(mrna, TARGET, 0) == 0
    assert distinct_matches(mrna, TARGET, 1) == 1


def test_two_mismatches_not_counted():
    two = "AGGUACGUACGUACGC"  # one substitution in each half
    mrna = _wrap(two)
    assert distinct_matches(mrna, TARGET, 1) == 0


def test_tandem_copies_counted_separately():
    # two full copies back to back: exactly one target length apart -> two clusters
    mrna = TARGET + TARGET
    assert distinct_matches(mrna, TARGET, 0) == 2


def test_overlapping_hits_collapse_to_one():
    # a target built from a period-8 repeat matches at 0 and 8; closer than one
    # target length apart, so they collapse into a single cluster
    target = "ACGUACGU" * 2  # 16 nt, self-similar with period 8
    mrna = "ACGUACGU" * 3  # matches at 0 and 8
    assert distinct_matches(mrna, target, 0) == 1


def test_compute_duplications_gates_and_raises():
    genes = np.array(["G1", "G2", "MISSING_LOCATED", "MISSING_UNLOCATED"])
    starts = np.array([0, -1, 5, -1])
    asos = np.array(["X"] * 4)  # unused for the missing-gene rows
    gene_mrna = {"G1": _wrap(TARGET), "G2": _wrap(TARGET)}

    # An unlocated ASO whose gene is also absent -> NaN, no error.
    de, dn = compute_duplications(starts[3:], genes[3:], asos[3:], gene_mrna)
    assert np.isnan(de[0]) and np.isnan(dn[0])

    # A located ASO whose gene sequence is missing -> loud failure.
    with pytest.raises(KeyError):
        compute_duplications(starts[2:3], genes[2:3], asos[2:3], gene_mrna)
