import os
import pytest
import pandas as pd
from unittest.mock import patch, MagicMock
from tauso.off_target.search import get_bowtie_index_base, run_bowtie_search, annotate_hits, find_all_gene_off_targets


# --- FIXTURES ---

@pytest.fixture
def mock_paths(tmp_path):
    """Mocks the data directory paths."""
    fasta_path = tmp_path / "GRCh38.fa"
    fasta_path.touch()  # Create dummy file

    with patch("tauso.off_target.search.get_paths") as mock:
        mock.return_value = {"fasta": str(fasta_path)}
        yield mock


@pytest.fixture
def mock_shutil_which():
    """Mocks shutil.which to pretend bowtie is installed."""
    with patch("shutil.which") as mock:
        mock.return_value = "/usr/bin/bowtie"
        yield mock


@pytest.fixture
def sample_sam_output():
    """
    Returns a sample SAM output string from Bowtie.
    Columns: QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, MRNM, MPOS, ISIZE, SEQ, QUAL, TAGS...
    """
    # 1. Perfect match on chr1, pos 100, + strand (Flag 0)
    line1 = "ASO_Seq\t0\tchr1\t100\t255\t20M\t*\t0\t0\tACGT\tIIII\tNM:i:0"
    # 2. 2 Mismatches on chr2, pos 500, - strand (Flag 16)
    line2 = "ASO_Seq\t16\tchr2\t500\t255\t20M\t*\t0\t0\tACGT\tIIII\tNM:i:2"
    # 3. Unmapped read (Flag 4) - Should be ignored
    line3 = "ASO_Seq\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII"

    return f"{line1}\n{line2}\n{line3}\n"


# --- TESTS ---

def test_get_bowtie_index_builds_if_missing(mock_paths, mock_shutil_which):
    """Test that bowtie-build is called if index doesn't exist."""
    with patch("os.path.exists", return_value=False), \
            patch("subprocess.run") as mock_run:
        index_base = get_bowtie_index_base()

        # Verify bowtie-build was called
        assert "bowtie-build" in mock_run.call_args[0][0]
        assert "GRCh38_bowtie_index" in index_base


def test_get_bowtie_index_skips_if_exists(mock_paths, mock_shutil_which):
    """Test that we skip building if index exists."""
    with patch("os.path.exists", return_value=True), \
            patch("subprocess.run") as mock_run:
        get_bowtie_index_base()

        # Should NOT call subprocess
        mock_run.assert_not_called()


def test_run_bowtie_search_parsing(mock_paths, mock_shutil_which, sample_sam_output):
    """Test parsing of SAM output into list of dicts."""
    with patch("tauso.off_target.search.get_bowtie_index_base", return_value="dummy_index"), \
            patch("subprocess.run") as mock_run:
        # Mock the subprocess stdout
        mock_process = MagicMock()
        mock_process.stdout = sample_sam_output
        mock_run.return_value = mock_process

        hits = run_bowtie_search("ACGT", max_mismatches=2)

        assert len(hits) == 2  # The unmapped read should be skipped

        # Check Hit 1 (+ strand, 0 mismatches)
        assert hits[0]['chrom'] == "chr1"
        assert hits[0]['start'] == 99  # 1-based SAM 100 -> 0-based 99
        assert hits[0]['strand'] == "+"
        assert hits[0]['mismatches'] == 0

        # Check Hit 2 (- strand, 2 mismatches)
        assert hits[1]['chrom'] == "chr2"
        assert hits[1]['start'] == 499
        assert hits[1]['strand'] == "-"
        assert hits[1]['mismatches'] == 2


def test_annotate_hits_priority():
    """Test that hits are annotated correctly with Exon > Intron > Gene priority."""

    # Mock Hits
    raw_hits = [
        {'chrom': 'chr1', 'start': 100, 'end': 120, 'mismatches': 0},
        {'chrom': 'chr2', 'start': 200, 'end': 220, 'mismatches': 1},
    ]

    # Mock Features to return from DB
    # Hit 1 lands in an EXON of Gene A
    feat_exon = MagicMock(featuretype='exon')
    feat_exon.attributes = {'gene_name': ['GeneA'], 'gene_id': ['ENSG01']}

    # Hit 1 ALSO lands in the GENE A (genomic span)
    feat_gene_a = MagicMock(featuretype='gene')
    feat_gene_a.attributes = {'gene_name': ['GeneA'], 'gene_id': ['ENSG01']}

    # Hit 2 lands ONLY in the GENE B (Intronic region, but no explicit intron feature found)
    # Or let's say it lands in an explicit Intron
    feat_intron = MagicMock(featuretype='intron')
    feat_intron.attributes = {'gene_name': ['GeneB'], 'gene_id': ['ENSG02']}

    # Mock the DB.region() call
    with patch("tauso.off_target.search.load_db") as mock_load_db:
        mock_db = MagicMock()
        # side_effect allows us to return different lists for different calls
        mock_db.region.side_effect = [
            [feat_gene_a, feat_exon],  # Result for Hit 1 (Order shouldn't matter)
            [feat_intron]  # Result for Hit 2
        ]
        mock_load_db.return_value = mock_db

        df = annotate_hits(raw_hits)

        assert len(df) == 2

        # Hit 1 should be 'exon' (Highest priority)
        assert df.iloc[0]['gene_name'] == 'GeneA'
        assert df.iloc[0]['region_type'] == 'exon'

        # Hit 2 should be 'intron'
        assert df.iloc[1]['gene_name'] == 'GeneB'
        assert df.iloc[1]['region_type'] == 'intron'


def test_annotate_hits_intergenic():
    """Test that hits with no DB overlap are labeled Intergenic."""
    raw_hits = [{'chrom': 'chr1', 'start': 1000, 'end': 1020}]

    with patch("tauso.data.load_db") as mock_load_db:
        mock_db = MagicMock()
        mock_db.region.return_value = []  # No features found
        mock_load_db.return_value = mock_db

        df = annotate_hits(raw_hits)

        assert df.iloc[0]['region_type'] == 'Intergenic'
        assert df.iloc[0]['gene_name'] is None


def test_full_pipeline_integration(mock_paths, mock_shutil_which, sample_sam_output):
    """Tests the find_all_gene_off_targets wrapper."""
    with patch("tauso.off_target.search.run_bowtie_search") as mock_search, \
            patch("tauso.off_target.search.annotate_hits") as mock_annotate:
        # Mock returns
        mock_search.return_value = [{'chrom': 'chr1', 'start': 100}]
        mock_annotate.return_value = pd.DataFrame([
            {'chrom': 'chr1', 'start': 100, 'gene_name': 'GeneA'}
        ])

        df = find_all_gene_off_targets("ACGT", max_mismatches=2)

        assert not df.empty
        assert 'gene_name' in df.columns
        mock_search.assert_called_once_with("ACGT", 2)
        mock_annotate.assert_called_once()