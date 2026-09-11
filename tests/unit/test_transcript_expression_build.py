"""What build-cohort-transcript-expression writes, given a cohort and the DepMap table.

The guards around this build are covered in test_transcript_expression_guards; here it is the
extraction itself: which row each cell line gets, which columns come out, and how the values
are derived.
"""

import json
import os

import pandas as pd
import pytest

from tauso.cli import build_cohort_transcript_expression

COHORT = {"CellA": "ACH-000001", "CellB": "ACH-000002", "Absent": "ACH-000999"}
TRANSCRIPTS = ["ENST00000000001.5", "ENST00000000002.1", "ENST00000000003.2"]


class _Feature:
    def __init__(self, accession, gene, name):
        self.id = accession
        self.attributes = {"gene_name": [gene], "transcript_name": [name]}


class _Db:
    def features_of_type(self, _types):
        return [
            _Feature("ENST00000000001.5", "GENEA", "GENEA-201"),
            _Feature("ENST00000000002.1", "GENEB", "GENEB-201"),
        ]


@pytest.fixture
def built(tmp_path, monkeypatch):
    """Run the build over a small table and return the output directory."""
    from tauso import cli

    rows = [
        # ACH-000001 has two profiles; only the default one may be used.
        {
            "ModelID": "ACH-000001",
            "IsDefaultEntryForModel": "No",
            TRANSCRIPTS[0]: 9.0,
            TRANSCRIPTS[1]: 9.0,
            TRANSCRIPTS[2]: 9.0,
        },
        {
            "ModelID": "ACH-000001",
            "IsDefaultEntryForModel": "Yes",
            TRANSCRIPTS[0]: 1.0,
            TRANSCRIPTS[1]: 3.0,
            TRANSCRIPTS[2]: 2.0,
        },
        {
            "ModelID": "ACH-000002",
            "IsDefaultEntryForModel": "Yes",
            TRANSCRIPTS[0]: 0.0,
            TRANSCRIPTS[1]: 0.0,
            TRANSCRIPTS[2]: 5.0,
        },
        {
            "ModelID": "ACH-000003",
            "IsDefaultEntryForModel": "Yes",
            TRANSCRIPTS[0]: 7.0,
            TRANSCRIPTS[1]: 7.0,
            TRANSCRIPTS[2]: 7.0,
        },
    ]
    pd.DataFrame(rows).to_csv(tmp_path / cli.TRANSCRIPT_EXPRESSION_CSV, index=False)
    (tmp_path / "cell_cohort.json").write_text(json.dumps(COHORT))

    monkeypatch.setattr(cli, "get_data_dir", lambda: str(tmp_path))
    monkeypatch.setattr(cli, "file_matches_hash", lambda *a, **k: True)
    monkeypatch.setattr(cli, "load_gtf_db", lambda *a, **k: _Db())
    build_cohort_transcript_expression.callback(force=True)
    return tmp_path / "processed_transcript_expression"


def test_one_file_per_cohort_cell_line_present_in_the_table(built):
    assert sorted(os.listdir(built)) == [
        ".cohort.json",
        "ACH-000001_transcript_expression.csv",
        "ACH-000002_transcript_expression.csv",
    ]


def test_a_cell_line_outside_the_cohort_is_not_written(built):
    assert not (built / "ACH-000003_transcript_expression.csv").exists()


def test_the_columns_are_fixed(built):
    d = pd.read_csv(built / "ACH-000001_transcript_expression.csv")
    assert list(d.columns) == ["Transcript", "TranscriptName", "Gene", "expression_norm", "expression_TPM"]


def test_the_default_profile_is_the_one_used(built):
    # ACH-000001's non-default row holds 9.0 everywhere; the default row holds 1, 3, 2.
    d = pd.read_csv(built / "ACH-000001_transcript_expression.csv")
    assert sorted(d.expression_norm) == [1.0, 2.0, 3.0]


def test_transcript_ids_lose_their_version(built):
    d = pd.read_csv(built / "ACH-000001_transcript_expression.csv")
    assert set(d.Transcript) == {"ENST00000000001", "ENST00000000002", "ENST00000000003"}


def test_tpm_is_the_log_undone(built):
    d = pd.read_csv(built / "ACH-000001_transcript_expression.csv")
    assert d.expression_TPM.tolist() == pytest.approx((2**d.expression_norm - 1).tolist())


def test_rows_come_out_most_expressed_first(built):
    d = pd.read_csv(built / "ACH-000002_transcript_expression.csv")
    assert d.expression_norm.tolist() == sorted(d.expression_norm, reverse=True)


def test_the_gene_and_name_come_from_the_annotation(built):
    d = pd.read_csv(built / "ACH-000001_transcript_expression.csv").set_index("Transcript")
    assert d.loc["ENST00000000001", "Gene"] == "GENEA"
    assert d.loc["ENST00000000001", "TranscriptName"] == "GENEA-201"
    # A transcript the annotation does not know keeps its row with no gene.
    assert pd.isna(d.loc["ENST00000000003", "Gene"])


def test_the_sentinel_names_the_cohort_built_for(built):
    assert json.loads((built / ".cohort.json").read_text()) == sorted(COHORT.values())
