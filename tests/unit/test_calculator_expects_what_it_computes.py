"""Every column a populate function fills must be on its step's expected list.

`_step` only saves the names it is told to expect, so a column that is computed but left
off the list is calculated on every run and then discarded -- silently, with no error and
no feature. That is what happened to the branch-point family.
"""

import re
from pathlib import Path

import pytest

CALCULATOR = Path(__file__).resolve().parents[2] / "src/tauso/populate/calculators/calculator.py"
POPULATE_STRUCTURE = Path(__file__).resolve().parents[2] / "src/tauso/populate/populate_structure.py"


def _expected_names_of_step(step):
    source = CALCULATOR.read_text()
    block = source[source.index(f"def {step}"):]
    block = block[: block.index("self._step")]
    return set(re.findall(r"STRUCTURE_SENSE_[A-Z_]+|STRUCT_SENSE_[A-Z_]+", block))


def _assigned_names_of_populate():
    source = POPULATE_STRUCTURE.read_text()
    return set(re.findall(r"all_data\[(STRUCTURE_SENSE_[A-Z_]+|STRUCT_SENSE_[A-Z_]+)\]", source))


def test_structure_step_expects_every_column_populate_assigns():
    assigned = _assigned_names_of_populate()
    expected = _expected_names_of_step("calculate_structure")
    assert assigned, "found no assigned columns; the parser needs updating"
    dropped = sorted(assigned - expected)
    assert not dropped, f"computed but never saved: {dropped}"


@pytest.mark.parametrize("name", sorted(_assigned_names_of_populate() & {
    "STRUCTURE_SENSE_BRANCH_POINT_OFFSET",
    "STRUCTURE_SENSE_BRANCH_POINT_SCORE",
    "STRUCTURE_SENSE_BRANCH_POINT_STRONG_COUNT",
    "STRUCTURE_SENSE_BRANCH_POINT_SIGNED_DIST",
    "STRUCTURE_SENSE_BRANCH_POINT_SIGNED_LOGDIST",
}))
def test_branch_point_columns_reach_the_step(name):
    assert name in _expected_names_of_step("calculate_structure")
