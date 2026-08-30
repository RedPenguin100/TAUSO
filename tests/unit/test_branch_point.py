import numpy as np
import pytest

from tauso.populate.populate_structure import _BRANCH_POINT_WINDOW, _find_branch_point

LO, HI = _BRANCH_POINT_WINDOW


def intron_with_branch_a(offset, before_two, before_one, after_one, length=200):
    """An intron whose only adenosine in the branch-point window sits `offset` nt from its end."""
    seq = list("GC" * (length // 2))
    # Blank the window first: its other positions are neighbours of the branch A, so filling them
    # afterwards would overwrite the motif being planted.
    for d in range(LO, HI + 1):
        seq[length - d] = "G"
    i = length - offset
    seq[i] = "A"
    seq[i - 2] = before_two
    seq[i - 1] = before_one
    seq[i + 1] = after_one
    return "".join(seq)


def test_finds_the_full_consensus():
    offset, score, strong = _find_branch_point(intron_with_branch_a(25, "C", "T", "C"))
    assert offset == 25.0
    assert score == 4.0
    assert strong == 1.0


def test_rna_and_dna_score_alike():
    dna = _find_branch_point(intron_with_branch_a(25, "C", "T", "C"))
    rna = _find_branch_point(intron_with_branch_a(25, "C", "U", "C").replace("T", "U"))
    assert dna == rna


def test_partial_consensus_scores_lower():
    _, score, strong = _find_branch_point(intron_with_branch_a(25, "G", "G", "G"))
    assert score == 1.0
    assert strong == 0.0


def test_prefers_the_better_motif_over_the_nearer_one():
    seq = list(intron_with_branch_a(40, "C", "T", "C"))
    seq[len(seq) - 20] = "A"
    offset, score, _ = _find_branch_point("".join(seq))
    assert offset == 40.0
    assert score == 4.0


@pytest.mark.parametrize("offset", [LO, HI])
def test_window_edges_are_included(offset):
    found, _, _ = _find_branch_point(intron_with_branch_a(offset, "C", "T", "C"))
    assert found == float(offset)


@pytest.mark.parametrize("offset", [LO - 1, HI + 1])
def test_outside_the_window_is_ignored(offset):
    seq = list("GC" * 100)
    for d in range(LO, HI + 2):
        seq[len(seq) - d] = "G"
    seq[len(seq) - offset] = "A"
    found, score, strong = _find_branch_point("".join(seq))
    assert np.isnan(found)
    assert np.isnan(score)
    assert strong == 0.0


def test_short_intron_has_no_branch_point():
    assert all(np.isnan(v) for v in _find_branch_point("ACGT" * 5))
    assert all(np.isnan(v) for v in _find_branch_point(None))


def test_counts_every_strong_candidate():
    seq = list("GC" * 100)
    for d in range(LO, HI + 1):
        seq[len(seq) - d] = "G"
    for d in (20, 30, 40):
        i = len(seq) - d
        seq[i], seq[i - 2], seq[i - 1], seq[i + 1] = "A", "C", "T", "C"
    _, score, strong = _find_branch_point("".join(seq))
    assert score == 4.0
    assert strong == 3.0
