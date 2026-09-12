"""build-cohort-transcript-expression over the real DepMap table and the real annotation.

Nothing here is stubbed but the directory the build reads and writes, so the transcript-to-gene
mapping and the row selection are the ones that run in production. The values are checked against
the CSV the table is published as, read independently of the build.
"""

import json
import os

import pandas as pd
import pytest

from tauso import cli
from tauso.cli import build_cohort_transcript_expression
from tauso.data.data import get_data_dir

pytestmark = pytest.mark.integration

# Both carry a default and a non-default profile, so they pin which row is chosen.
COHORT = {"CellA": "ACH-000004", "CellB": "ACH-000008", "Absent": "ACH-000000"}


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    source = os.path.join(get_data_dir(), cli.TRANSCRIPT_EXPRESSION_PARQUET)
    data_dir = tmp_path_factory.mktemp("transcript_build")
    os.symlink(source, data_dir / cli.TRANSCRIPT_EXPRESSION_PARQUET)
    (data_dir / "cell_cohort.json").write_text(json.dumps(COHORT))

    original = cli.get_data_dir
    cli.get_data_dir = lambda: str(data_dir)
    try:
        build_cohort_transcript_expression.callback(force=True)
    finally:
        cli.get_data_dir = original
    return data_dir / "processed_transcript_expression"


@pytest.fixture(scope="module")
def source_rows():
    """The same two cell lines read straight from the table, as the answer to check against."""
    source = os.path.join(get_data_dir(), cli.TRANSCRIPT_EXPRESSION_CSV)
    ids = pd.read_csv(source, usecols=["ModelID"])["ModelID"]
    wanted = set(ids.index[ids.isin(COHORT.values())])
    d = pd.read_csv(source, skiprows=lambda i: i > 0 and (i - 1) not in wanted)
    return d[d["IsDefaultEntryForModel"] == "Yes"].set_index("ModelID")


def test_one_file_per_cohort_cell_line_in_the_table(built):
    assert sorted(os.listdir(built)) == [
        ".cohort.json",
        "ACH-000004_transcript_expression.csv",
        "ACH-000008_transcript_expression.csv",
    ]


@pytest.mark.parametrize("model_id", ["ACH-000004", "ACH-000008"])
def test_the_values_are_the_default_profile_row(built, source_rows, model_id):
    out = pd.read_csv(built / f"{model_id}_transcript_expression.csv").set_index("Transcript")
    want = source_rows.loc[model_id]
    columns = [c for c in want.index if c.startswith("ENST")]
    assert len(out) == len(columns)
    for column in columns[:200]:
        assert out.loc[column.split(".", 1)[0], "expression_norm"] == pytest.approx(float(want[column]))


def test_tpm_is_the_log_undone(built):
    d = pd.read_csv(built / "ACH-000004_transcript_expression.csv")
    assert d.expression_TPM.tolist() == pytest.approx((2**d.expression_norm - 1).tolist())


def test_rows_come_out_most_expressed_first(built):
    d = pd.read_csv(built / "ACH-000004_transcript_expression.csv")
    assert d.expression_norm.tolist() == sorted(d.expression_norm, reverse=True)


def test_the_gene_names_come_from_the_annotation(built):
    from tauso.data.data import load_gtf_db

    d = pd.read_csv(built / "ACH-000004_transcript_expression.csv").dropna(subset=["Gene"])
    named = d.set_index("Transcript")["Gene"].to_dict()
    checked = 0
    for feature in load_gtf_db().features_of_type(("mRNA", "transcript")):
        accession = feature.id.split(".", 1)[0]
        if accession in named and feature.attributes.get("gene_name"):
            assert named[accession] == feature.attributes["gene_name"][0]
            checked += 1
            if checked == 200:
                return
    assert checked, "no transcript in the output matched the annotation"


def test_the_columns_are_fixed(built):
    d = pd.read_csv(built / "ACH-000004_transcript_expression.csv", nrows=1)
    assert list(d.columns) == ["Transcript", "TranscriptName", "Gene", "expression_norm", "expression_TPM"]


def test_the_sentinel_names_the_cohort_built_for(built):
    assert json.loads((built / ".cohort.json").read_text()) == sorted(COHORT.values())
