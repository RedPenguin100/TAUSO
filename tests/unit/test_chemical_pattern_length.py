"""A chemical pattern that does not have one character per residue is a corrupt row.

The pattern is positional, so a length that disagrees with the sequence means every sugar read
out of it may belong to a different residue. Scoring NaN hides that; the features raise instead.
"""

import numpy as np
import pytest

from tauso.common.modifications import check_pattern_length
from tauso.features.hybridization.hybridization_features import get_dna_rna_dg_region
from tauso.features.hybridization.md_weights import get_moe_md_contribution

MOE_GAPMER = "MMMMMddddddddddMMMMM"
TWENTY_MER = "ACGTACGTACGTACGTACGT"
MOE_LABEL = "MOE/5-methylcytosines/deoxy"


def test_a_matching_pattern_is_accepted():
    check_pattern_length(TWENTY_MER, MOE_GAPMER)


@pytest.mark.parametrize("sequence", [TWENTY_MER[:18], TWENTY_MER + "AC"])
def test_a_pattern_of_the_wrong_length_raises(sequence):
    with pytest.raises(ValueError, match="chemical_pattern length"):
        check_pattern_length(sequence, MOE_GAPMER)


def test_a_missing_pattern_is_left_to_the_caller():
    """Absent is not the same as wrong: the caller decides what an oligo with no pattern means."""
    check_pattern_length(TWENTY_MER, None)
    check_pattern_length(TWENTY_MER, np.nan)


@pytest.mark.parametrize(
    ("name", "call"),
    [
        ("hybridization", lambda seq: get_dna_rna_dg_region(seq, MOE_GAPMER, "wing5")),
        ("moe md", lambda seq: get_moe_md_contribution(seq, MOE_GAPMER, MOE_LABEL)),
    ],
)
def test_the_hybridization_features_raise_rather_than_score_nan(name, call):
    assert np.isfinite(call(TWENTY_MER)), name
    with pytest.raises(ValueError, match="chemical_pattern length"):
        call(TWENTY_MER[:18])


@pytest.mark.parametrize(
    ("name", "call"),
    [
        ("hybridization", lambda p: get_dna_rna_dg_region(TWENTY_MER, p, "wing5")),
        ("moe md", lambda p: get_moe_md_contribution(TWENTY_MER, p, MOE_LABEL)),
    ],
)
def test_no_pattern_still_scores_nan(name, call):
    assert np.isnan(call(None)), name
