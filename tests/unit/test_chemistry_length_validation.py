"""The chemistry strings must line up with the sequence.

`chemical_pattern` is one character per residue, `ps_pattern` one per linkage. Both hold
by construction, so a violation means a parser met a modification it did not know and
dropped it -- the failure mode behind the 651 mis-parsed 2'-OMe rows.
"""

import pandas as pd
import pytest
from notebooks.models.utility import validate_chemistry_lengths

from tauso.data.consts import ASO_SEQUENCE, CHEMICAL_PATTERN, PS_PATTERN

INDEX = "index_oligo"


def frame(rows):
    return pd.DataFrame(rows, columns=[INDEX, ASO_SEQUENCE, CHEMICAL_PATTERN, PS_PATTERN])


GOOD = [
    (1, "ACGTACGT", "MMddddMM", "sssssss"),
    (2, "ACGT", "CCdd", "sss"),
]


def test_well_formed_rows_pass():
    validate_chemistry_lengths(frame(GOOD), INDEX)


def test_short_chemical_pattern_raises_and_names_the_row():
    rows = GOOD + [(77, "ACGTACGT", "MMdddMM", "sssssss")]
    with pytest.raises(ValueError) as e:
        validate_chemistry_lengths(frame(rows), INDEX)
    assert CHEMICAL_PATTERN in str(e.value)
    assert "77" in str(e.value)


def test_ps_pattern_must_be_one_shorter_than_the_sequence():
    rows = GOOD + [(88, "ACGT", "CCdd", "ssss")]
    with pytest.raises(ValueError) as e:
        validate_chemistry_lengths(frame(rows), INDEX)
    assert PS_PATTERN in str(e.value)
    assert "88" in str(e.value)


def test_a_pattern_equal_in_length_to_ps_pattern_is_caught():
    """The 651-row signature: a chemical_pattern one short reads as a valid ps_pattern."""
    rows = [(99, "ACGTACGTACGTACGT", "CCCddddddddddCC", "sssssssssssssss")]
    with pytest.raises(ValueError) as e:
        validate_chemistry_lengths(frame(rows), INDEX)
    assert "99" in str(e.value)


def test_missing_pattern_is_a_violation():
    rows = GOOD + [(55, "ACGT", None, "sss")]
    with pytest.raises(ValueError):
        validate_chemistry_lengths(frame(rows), INDEX)


def test_the_message_caps_the_list_but_reports_the_true_count():
    rows = [(i, "ACGT", "CC", "sss") for i in range(15)]
    with pytest.raises(ValueError) as e:
        validate_chemistry_lengths(frame(rows), INDEX)
    assert "15 of 15" in str(e.value)
    assert "+5 more" in str(e.value)
