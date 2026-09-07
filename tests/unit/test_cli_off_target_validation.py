"""`run-off-target` must reject a sequence bowtie would happily mis-search.

Validation runs before the genome check, so these need no genome installed. An accepted
sequence is covered by test_mipomersen_real_search, which runs the real search.
"""

import pytest
from click.testing import CliRunner

from tauso.cli import main


def run(sequence):
    return CliRunner().invoke(main, ["run-off-target", sequence])


@pytest.mark.parametrize(
    "sequence, expected",
    [
        ("ACGTNACGT", "non-ACGU"),
        ("ACGT-ACGT", "non-ACGU"),
        ("ACGTRACGT", "non-ACGU"),  # IUPAC ambiguity code
        ("ACGT ACGT", "whitespace"),
        ("", "empty"),
    ],
)
def test_bad_sequences_are_rejected(sequence, expected):
    result = run(sequence)
    assert result.exit_code == 1
    assert "Invalid sequence" in result.output
    assert expected in result.output
