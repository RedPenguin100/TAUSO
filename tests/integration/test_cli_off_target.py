import pandas as pd
import pytest

from click.testing import CliRunner
from tauso.cli import main



@pytest.mark.integration
def test_mipomersen_real_search():
    """
    Integration Test: Runs Mipomersen (RNA) against the real local genome.

    Requires: 'tauso setup-genome' to have been run previously.
    Verifies:
        1. RNA (U) -> DNA (T) normalization occurs.
        2. The biological target (APOB on chr2) is found.
    """
    runner = CliRunner()

    # Mipomersen sequence with Uracils
    mipomersen_rna = "GCCUCAGTCTGCTTCGCACC"

    # Run the real command (no mocks)
    result = runner.invoke(main, ["run-off-target", "-o", "off_target.csv", mipomersen_rna])
    df = pd.read_csv('off_target.csv')


    # If the genome isn't installed, this might fail with exit_code 1
    if result.exit_code != 0:
        pytest.fail(f"CLI failed. Output:\n{result.output}")

    assert result.exit_code == 0


    assert "APOB" in result.output
    assert "chr2" in result.output
    assert "exon" in result.output

    assert (df["mismatches"] == 0).all()
