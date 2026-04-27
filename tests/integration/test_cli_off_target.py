import pandas as pd
import pytest
from click.testing import CliRunner

from tauso.cli import main
from tauso.off_target.search import (
    annotate_hits,
    run_bowtie_search,
)


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
    df = pd.read_csv("off_target.csv")

    # If the genome isn't installed, this might fail with exit_code 1
    if result.exit_code != 0:
        pytest.fail(f"CLI failed. Output:\n{result.output}")

    assert result.exit_code == 0

    assert "APOB" in result.output
    assert "chr2" in result.output
    assert "exon" in result.output

    # Changed from .all() to .any() to verify the perfect match exists
    # among the off-targets.
    assert (df["mismatches"] == 0).any()


def test_inotersen_off_target_search():
    """
    Real end-to-end integration test for off-target search.
    Uses Inotersen (Tegsedi), a 20-mer ASO that targets the TTR gene.
    """
    # The actual sequence of Inotersen
    inotersen_seq = "TCTTGGTTACATGAAATCCC"

    # 1. Run the real Bowtie search against the local GRCh38 index
    # We allow up to 2 mismatches to see the primary target and potential off-targets
    hits, total_counts = run_bowtie_search(inotersen_seq, max_mismatches=2)

    # Ensure Bowtie actually executed and found the primary target
    assert len(hits) > 0, "Bowtie failed to find any hits in the genome."

    # 2. Run the real DB annotation (Queries the gffutils/sqlite database)
    annotated_df = annotate_hits(hits)

    # Verify the output is a properly populated DataFrame
    assert isinstance(annotated_df, pd.DataFrame)
    assert not annotated_df.empty

    # 3. Validate the biological correctness of the result
    genes_found = annotated_df["gene_name"].tolist()

    # Inotersen MUST hit the TTR gene
    assert "TTR" in genes_found, f"Expected to find TTR gene, but got: {genes_found}"

    # 4. Check the accuracy of the primary target
    ttr_hits = annotated_df[annotated_df["gene_name"] == "TTR"]

    # The primary target should be a perfect match (0 mismatches)
    assert ttr_hits["mismatches"].min() == 0, "Expected a perfect match for the TTR gene."

    # Ensure the region logic works (it should map to an exon or intron of TTR)
    assert ttr_hits["region_type"].iloc[0] in ["exon", "intron", "transcript", "gene"]
