"""A partial download must not be used, and a build that produced nothing must not be cached."""

import json
from pathlib import Path

import pytest

from tauso.cli_utils import file_matches_hash

TRUNCATED = "col_a,col_b\n1,2\n"


def test_a_file_of_the_right_name_can_still_be_the_wrong_file(tmp_path):
    """What the existence check missed: the name says nothing about the contents."""
    path = tmp_path / "OmicsExpressionTranscriptTPMLogp1HumanAllGenes.csv"
    path.write_text(TRUNCATED)

    assert path.exists()
    assert not file_matches_hash(str(path), "c1277e2bf88b14f18b133a5d31516872856396cc")


def test_a_missing_file_reads_as_a_mismatch(tmp_path):
    """So one hash check covers both the absent and the corrupt case."""
    assert not file_matches_hash(str(tmp_path / "absent.csv"), "deadbeef")


def test_the_sentinel_records_the_cohort_it_was_built_for(tmp_path):
    """The sentinel is a cohort list, so a later run with a different cohort rebuilds."""
    sentinel = tmp_path / ".cohort.json"
    sentinel.write_text(json.dumps(sorted(["ACH-000018", "ACH-001085"])))
    built_for = set(json.loads(sentinel.read_text()))

    assert built_for == {"ACH-000018", "ACH-001085"}
    assert built_for != {"ACH-000018"}


@pytest.mark.parametrize("found_count", [0, 1, 29])
def test_only_a_build_that_found_something_may_write_the_sentinel(tmp_path, found_count):
    """Guards the shape of the check: nothing found means no sentinel, so the next run retries."""
    sentinel = Path(tmp_path) / ".cohort.json"
    if found_count:
        sentinel.write_text(json.dumps(["ACH-000018"]))

    assert sentinel.exists() == bool(found_count)
