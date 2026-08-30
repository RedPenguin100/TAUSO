"""`window_starts` decides which windows of the fold each anchoring averages over.

The orientation is the part worth pinning: the oligo binds antiparallel, so its 5'
terminus pairs with the target's 3' side. Getting that backwards silently swaps two
features rather than raising, and the frozen regressions would only say "numbers moved".
"""

import pytest

from tauso.features.fold.vienna_access import window_starts

SENSE_START, SENSE_LENGTH = 60, 20


def starts(anchor, open_len):
    return list(window_starts(anchor, SENSE_START, SENSE_LENGTH, open_len))


@pytest.mark.parametrize("open_len", [4, 6, 13, 32])
def test_sweeps_give_one_window_per_target_position(open_len):
    assert len(starts("a5", open_len)) == SENSE_LENGTH
    assert len(starts("a3", open_len)) == SENSE_LENGTH


@pytest.mark.parametrize("open_len", [4, 6, 13, 32])
def test_a5_starts_inside_the_target_and_a3_ends_inside_it(open_len):
    a5, a3 = starts("a5", open_len), starts("a3", open_len)
    assert a5[0] == SENSE_START
    assert a5[-1] == SENSE_START + SENSE_LENGTH - 1
    assert a3[0] + open_len - 1 == SENSE_START
    assert a3[-1] + open_len - 1 == SENSE_START + SENSE_LENGTH - 1


@pytest.mark.parametrize("open_len", [4, 6, 13])
def test_ends_are_single_windows_flush_with_each_target_edge(open_len):
    (three,) = starts("aso3end", open_len)
    (five,) = starts("aso5end", open_len)
    # ASO 3' pairs with the target 5' side, so its window opens at the target's start.
    assert three == SENSE_START
    # ASO 5' pairs with the target 3' side, so its window closes at the target's end.
    assert five + open_len == SENSE_START + SENSE_LENGTH


def test_end_windows_are_the_extremes_of_the_sweeps():
    """Each end is a window the sweeps already visit, which is why they cost nothing."""
    open_len = 6
    assert starts("aso3end", open_len)[0] == starts("a5", open_len)[0]
    assert starts("aso5end", open_len)[0] == starts("a3", open_len)[-1]


def test_unknown_anchor_raises():
    with pytest.raises(ValueError, match="unknown anchor"):
        window_starts("middle", SENSE_START, SENSE_LENGTH, 6)


def test_reducers_are_named_and_unknown_ones_raise():
    from tauso.features.fold.vienna_access import reduce_energies

    energies = [1.0, 2.0, 6.0]
    assert reduce_energies(energies, "mean") == 3.0
    assert reduce_energies(energies, "std") == pytest.approx(2.160246899469287)
    with pytest.raises(ValueError, match="unknown reducer"):
        reduce_energies(energies, "median")


def test_a_single_window_has_no_spread():
    """Why the end anchorings are mean-only: std over one window is always zero."""
    from tauso.features.fold.vienna_access import reduce_energies

    assert reduce_energies([4.2], "std") == 0.0


def test_feature_name_marks_std_and_leaves_mean_unmarked():
    from tauso.populate.populate_fold import access_feature_name

    assert access_feature_name(60, None, 6, "a5", "mean") == "access_f60_sinf_u6_a5"
    assert access_feature_name(60, None, 6, "a5", "std") == "access_f60_sinf_u6_a5_std"
    with pytest.raises(ValueError, match="unknown reducer"):
        access_feature_name(60, None, 6, "a5", "median")
