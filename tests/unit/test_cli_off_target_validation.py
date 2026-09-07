"""`run-off-target` must reject a sequence bowtie would happily mis-search.

Validation runs before the genome check, so a rejected sequence needs no genome. An
accepted one goes on to the search, which is a transcriptome-wide bowtie run holding the
genome in memory, so the search is stubbed out: what is under test is the validation.
"""

import pytest
from click.testing import CliRunner

import tauso.cli as cli
from tauso.cli import main


def run(sequence):
    return CliRunner().invoke(main, ["run-off-target", sequence])


@pytest.fixture
def no_search(monkeypatch):
    """Let an accepted sequence reach the search without running it."""
    import pandas as pd

    monkeypatch.setattr(cli, "find_all_gene_off_targets", lambda *a, **k: pd.DataFrame())
    return monkeypatch


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


def test_surrounding_whitespace_is_stripped_not_rejected(no_search):
    """A quoted shell argument with a stray space is a paste artifact, not a bad sequence."""
    result = run("  ACGTACGT  ")
    assert "Invalid sequence" not in result.output


def test_rna_is_still_accepted_and_announced(no_search):
    result = run("ACGUACGU")
    assert "Invalid sequence" not in result.output
    assert "Normalized input sequence: ACGUACGU -> ACGTACGT" in result.output
