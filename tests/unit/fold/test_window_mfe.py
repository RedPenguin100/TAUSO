"""`window_mfe` reads a window's MFE out of a fold of some enclosing sequence.

Every test here pins the one invariant the whole shared-fold scheme rests on:

    window_mfe(sequence, start, size, span) == RNA.fold(sequence[start:start + size])[1]

with exact float equality, not a tolerance. The frozen feature regressions would also
catch a break, but only as "1000 numbers moved"; these say which window and by how much.
"""

import random

import numpy as np
import pytest
import ViennaRNA as RNA

from tauso.features.fold._vienna_internal import window_mfe
from tauso.features.fold.vienna_fold import (
    SharedFold,
    calculate_avg_mfe,
    calculate_avg_mfe_per_setting,
    calculate_end_mfe,
)

# The window sizes populate_fold actually asks for (DEFAULT_SETTINGS).
PRODUCTION_WINDOWS = [25, 40, 55, 70, 100, 150]


def random_rna(length, seed, alphabet="ACGU"):
    rng = random.Random(seed)
    return "".join(rng.choice(alphabet) for _ in range(length))


def own_region(mrna, start, sense_length, flank):
    """The cut a setting would fold on its own, with nothing shared from neighbouring sites."""
    return max(0, start - flank), min(len(mrna), start + sense_length + flank)


def assert_matches_direct_fold(sequence, start, window_size, max_bp_span=None):
    """window_mfe agrees exactly with folding that window on its own.

    `max_bp_span` of None means the tightest span that is still valid, the window itself.
    """
    expected = RNA.fold(sequence[start : start + window_size])[1]
    actual = window_mfe(sequence, start, window_size, max_bp_span or window_size)
    assert actual == expected, (
        f"window [{start}:{start + window_size}] of a {len(sequence)}nt sequence: "
        f"got {actual!r}, folding it alone gives {expected!r} (difference {actual - expected!r})"
    )


@pytest.mark.parametrize("window_size", PRODUCTION_WINDOWS)
@pytest.mark.parametrize("offset", [0, 1, 7, 40])
def test_matches_direct_fold_at_production_window_sizes(window_size, offset):
    sequence = random_rna(260, seed=window_size + offset)
    assert_matches_direct_fold(sequence, offset, window_size)


@pytest.mark.parametrize("seed", range(12))
def test_matches_direct_fold_on_random_sequences(seed):
    sequence = random_rna(180, seed=1000 + seed)
    rng = random.Random(seed)
    window_size = rng.choice(PRODUCTION_WINDOWS[:-1])
    start = rng.randrange(0, len(sequence) - window_size + 1)
    assert_matches_direct_fold(sequence, start, window_size)


@pytest.mark.parametrize("window_size", [25, 70, 150])
def test_matches_at_both_sequence_edges(window_size):
    """A stem at a window edge has no dangling neighbour; the edges are where that shows."""
    sequence = random_rna(200, seed=window_size)
    assert_matches_direct_fold(sequence, 0, window_size)
    assert_matches_direct_fold(sequence, len(sequence) - window_size, window_size)


def test_matches_when_the_window_is_the_whole_sequence():
    sequence = random_rna(120, seed=5)
    assert_matches_direct_fold(sequence, 0, len(sequence))


@pytest.mark.parametrize("window_size", [5, 6, 7, 8, 12])
def test_matches_for_windows_barely_longer_than_a_hairpin(window_size):
    """Windows near the minimum loop size, where the recursion's bounds are tightest."""
    sequence = random_rna(60, seed=window_size)
    for start in (0, 3, 60 - window_size):
        assert_matches_direct_fold(sequence, start, window_size)


@pytest.mark.parametrize("alphabet", ["A", "AC", "GC", "AU"])
def test_matches_on_homopolymers_and_restricted_alphabets(alphabet):
    """Sequences where no pair, or only one kind of pair, can form."""
    sequence = random_rna(120, seed=len(alphabet), alphabet=alphabet)
    for window_size in (25, 70):
        assert_matches_direct_fold(sequence, 10, window_size)


def test_unpairable_sequence_has_zero_mfe():
    """Nothing can pair in a poly-A window, so the MFE is the empty structure's."""
    sequence = "A" * 100
    assert window_mfe(sequence, 0, 50, 50) == 0.0


@pytest.mark.parametrize("window_size", [25, 100])
def test_matches_on_sequences_containing_ambiguous_bases(window_size):
    """Genomic sequence carries Ns; they must encode to "pairs with nothing", as ViennaRNA does."""
    sequence = random_rna(200, seed=window_size, alphabet="ACGUACGUACGUN")
    for start in (0, 17, 200 - window_size):
        assert_matches_direct_fold(sequence, start, window_size)


@pytest.mark.parametrize("window_size", [25, 70])
@pytest.mark.parametrize("max_bp_span", [None, 150, 200])
def test_result_does_not_depend_on_the_span_cap(window_size, max_bp_span):
    """Any span at or above the window size is exact, because nothing inside a window
    can pair further apart than the window is wide."""
    sequence = random_rna(200, seed=window_size)
    assert_matches_direct_fold(sequence, 25, window_size, max_bp_span)


def test_every_production_window_size_out_of_one_shared_fold():
    """The whole point: one fold, every window size, all still exact."""
    sequence = random_rna(300, seed=42)
    widest = max(PRODUCTION_WINDOWS)
    for window_size in PRODUCTION_WINDOWS:
        for start in range(0, len(sequence) - window_size + 1, 37):
            assert_matches_direct_fold(sequence, start, window_size, widest)


@pytest.mark.parametrize("window_size,step", [(25, 4), (70, 5), (100, 7)])
def test_shared_fold_gives_the_same_answer_as_folding_the_cut_alone(window_size, step):
    """calculate_avg_mfe must not care which enclosing sequence the fold came from."""
    enclosing = random_rna(300, seed=window_size)
    offset, cut_length = 40, 220
    cut = enclosing[offset : offset + cut_length]

    standalone = calculate_avg_mfe(60, 20, window_size, step, fold=SharedFold(cut, 0, len(cut), window_size))
    shared = calculate_avg_mfe(
        60, 20, window_size, step, fold=SharedFold(enclosing, offset, len(cut), max(PRODUCTION_WINDOWS))
    )
    assert shared == standalone


def test_per_setting_matches_running_each_setting_on_its_own():
    """calculate_avg_mfe_per_setting is only a scheduling change over per-setting calls."""
    settings = [(30, 25, 4), (30, 40, 5), (60, 55, 5), (60, 70, 5), (120, 100, 7), (120, 150, 10)]
    mrna = random_rna(1200, seed=7)
    global_start, sense_length = 500, 20

    widest_flank = max(flank for flank, _, _ in settings)
    combined = calculate_avg_mfe_per_setting(
        mrna, global_start, sense_length, settings, own_region(mrna, global_start, sense_length, widest_flank)
    )

    for flank_size, window_size, step in settings:
        cut_start = max(0, global_start - flank_size)
        cut_end = min(len(mrna), global_start + sense_length + flank_size)
        cut = mrna[cut_start:cut_end]
        expected = calculate_avg_mfe(
            global_start - cut_start,
            sense_length,
            window_size,
            step,
            fold=SharedFold(cut, 0, len(cut), window_size),
        )
        assert combined[(flank_size, window_size, step)] == expected, (
            f"setting (flank={flank_size}, window={window_size}, step={step}) disagrees"
        )


def test_setting_is_nan_when_the_cut_is_shorter_than_its_window():
    """A window wider than the available sequence yields no value rather than raising."""
    settings = [(30, 150, 10)]
    mrna = random_rna(60, seed=3)
    result = calculate_avg_mfe_per_setting(mrna, 10, 20, settings, own_region(mrna, 10, 20, 30))
    assert np.isnan(result[(30, 150, 10)])


def test_a_wider_shared_region_does_not_change_any_setting():
    """populate_fold folds one region spanning several neighbouring target sites.

    Each setting must still be swept over its own cut. Bounding the sweep by the folded
    region instead would let the narrow-flank settings run past their cut into a
    neighbour's sequence, which no golden covers because they all fold the default region.
    """
    settings = [(30, 25, 4), (30, 40, 5), (60, 55, 5), (60, 70, 5), (120, 100, 7), (120, 150, 10)]
    mrna = random_rna(1200, seed=11)
    global_start, sense_length = 500, 20

    widest_flank = max(flank for flank, _, _ in settings)
    narrow = calculate_avg_mfe_per_setting(
        mrna, global_start, sense_length, settings, own_region(mrna, global_start, sense_length, widest_flank)
    )
    wide = calculate_avg_mfe_per_setting(mrna, global_start, sense_length, settings, fold_region=(200, 1000))
    assert narrow == wide


def test_a_wider_shared_region_does_not_change_the_end_features():
    mrna = random_rna(1200, seed=13)
    narrow = calculate_end_mfe(mrna, 500, 20, 30, 40, 5, 6, own_region(mrna, 500, 20, 30))
    wide = calculate_end_mfe(mrna, 500, 20, 30, 40, 5, 6, fold_region=(200, 1000))
    assert narrow == wide


@pytest.mark.parametrize(
    "offset,subseq_length",
    [(150, 120), (0, 201), (-1, 50)],
)
def test_a_subsequence_outside_the_fold_is_rejected(offset, subseq_length):
    """Reading past the folded region returns garbage energies rather than raising, so the
    sub-sequence bounds are checked where they are set instead of where they are used."""
    sequence = random_rna(200, seed=4)
    with pytest.raises(ValueError, match="does not fit"):
        SharedFold(sequence, offset, subseq_length, 70)


def test_a_subsequence_filling_the_whole_fold_is_allowed():
    sequence = random_rna(200, seed=4)
    assert SharedFold(sequence, 0, 200, 70).subseq_length == 200


@pytest.mark.parametrize("window_size,max_bp_span", [(150, 40), (70, 25), (25, 24)])
def test_a_window_wider_than_the_span_is_rejected(window_size, max_bp_span):
    """Pairs wider than the span were never computed, so the window would come back
    under-structured -- a plausible-looking energy several kcal/mol too high."""
    sequence = random_rna(400, seed=11)
    with pytest.raises(ValueError, match="max_bp_span"):
        window_mfe(sequence, 0, window_size, max_bp_span)
