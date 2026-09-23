"""A DepMap table written by depmap_table reads back as the CSV it came from, however it is stored."""

from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.csv as pacsv
import pytest

from tauso import cli
from tauso.cli_utils import sha256_file
from tauso.expression import depmap_table
from tauso.expression.depmap_table import (
    column_means,
    column_names,
    read_rows,
    record_sha256,
    recorded_sha256,
    table_exists,
    table_sha256,
    write_table,
)
from tauso.expression.general import get_general_expression_of_genes

HEADER = ",SequencingID,ModelID,IsDefaultEntryForModel,ModelConditionID,IsDefaultEntryForMC"
PROFILES = [
    ("0", "S1", "ACH-000001", "Yes", "MC1", "Yes"),
    ("1", "S2", "ACH-000002", "No", "MC2", "No"),
    ("2", "S3", "ACH-000002", "Yes", "MC3", "Yes"),
    ("3", "S4", "ACH-000003", "Yes", "MC4", "Yes"),
]
TRANSCRIPTS = [f"ENST{n:011d}.1" for n in range(7)]
GENES = ["TSPAN6 (7105)", "TNMD (64102)", "DPM1 (8813)", "SCYL3 (57147)", "FGR (2268)", "ENSG00000288724.1"]
TABLES = {"pyarrow": (cli.TRANSCRIPT_EXPRESSION_CSV, TRANSCRIPTS), "pandas": (cli.GENE_EXPRESSION_CSV, GENES)}


def write_csv(directory, name, columns):
    rng = np.random.default_rng(0)
    lines = [HEADER + "," + ",".join(columns)]
    for profile in PROFILES:
        values = [repr(float(v)) for v in rng.random(len(columns)) * 12]
        if profile[0] == "0":
            values[2] = ""  # a missing value
        lines.append(",".join(profile) + "," + ",".join(values))
    path = Path(directory) / name
    path.write_text("\n".join(lines) + "\n")
    return path


def parsed(csv, parser):
    """The whole CSV, parsed as the table's parser parses it."""
    if parser == "pandas":
        return pd.read_csv(csv)
    return pacsv.read_csv(csv).to_pandas()


@pytest.fixture(autouse=True)
def small_parts(monkeypatch):
    monkeypatch.setattr(depmap_table, "COLUMNS_PER_PART", 2)


@pytest.mark.parametrize("parser", TABLES)
def test_the_table_reads_back_as_the_csv(tmp_path, parser):
    name, columns = TABLES[parser]
    csv = write_csv(tmp_path, name, columns)
    whole = parsed(csv, parser)

    write_table(csv, parser)
    ids, got_columns, values = read_rows(csv, {"ACH-000001", "ACH-000002"})

    assert not csv.exists() and table_exists(csv)
    assert column_names(csv) == got_columns == columns
    assert list(ids) == ["ACH-000001", "ACH-000002"]
    # ACH-000002 is read from its default profile, the third row; the missing value reads as 0.
    want = whole.iloc[[0, 2]][columns].fillna(0.0).to_numpy(dtype=np.float64)
    np.testing.assert_array_equal(values, want)


@pytest.mark.parametrize("parser", TABLES)
def test_column_means_are_in_the_order_asked(tmp_path, parser):
    name, columns = TABLES[parser]
    csv = write_csv(tmp_path, name, columns)
    whole = parsed(csv, parser)
    write_table(csv, parser)

    asked = [columns[4], columns[0], columns[2]]
    np.testing.assert_array_equal(column_means(csv, asked), [np.nanmean(whole[c].to_numpy(float)) for c in asked])


def test_an_older_single_file_reads_the_same(tmp_path):
    (tmp_path / "parts").mkdir()
    (tmp_path / "single").mkdir()
    parts_csv = write_csv(tmp_path / "parts", cli.GENE_EXPRESSION_CSV, GENES)
    single_csv = write_csv(tmp_path / "single", cli.GENE_EXPRESSION_CSV, GENES)
    write_table(parts_csv, "pandas")
    pd.read_csv(single_csv).to_parquet(single_csv.with_suffix(".parquet"), index=False)

    ids = {"ACH-000001", "ACH-000002", "ACH-000003"}
    got, want = read_rows(parts_csv, ids), read_rows(single_csv, ids)
    assert list(got[0]) == list(want[0]) and got[1] == want[1]
    np.testing.assert_array_equal(got[2], want[2])
    np.testing.assert_array_equal(column_means(parts_csv, GENES), column_means(single_csv, GENES))


def test_writing_replaces_an_older_single_file(tmp_path):
    csv = write_csv(tmp_path, cli.GENE_EXPRESSION_CSV, GENES)
    single = csv.with_suffix(".parquet")
    pd.read_csv(csv).to_parquet(single, index=False)
    record_sha256(csv)

    write_table(csv, "pandas")

    assert not single.exists() and not Path(f"{single}.sha256").exists()
    assert recorded_sha256(csv) == table_sha256(csv)


def test_a_killed_conversion_is_not_reused(tmp_path):
    csv = write_csv(tmp_path, cli.TRANSCRIPT_EXPRESSION_CSV, TRANSCRIPTS)
    partial = tmp_path / (csv.stem + "_parts.partial")
    partial.mkdir()
    (partial / "part-000.parquet").write_text("half")
    assert not table_exists(csv)

    write_table(csv, "pyarrow")

    assert not partial.exists()
    assert read_rows(csv, {"ACH-000003"})[1] == TRANSCRIPTS


def test_the_recorded_hash_notices_a_changed_part(tmp_path):
    csv = write_csv(tmp_path, cli.GENE_EXPRESSION_CSV, GENES)
    write_table(csv, "pandas")
    assert recorded_sha256(csv) == table_sha256(csv)

    part = sorted((tmp_path / (csv.stem + "_parts")).glob("part-*.parquet"))[1]
    part.write_bytes(part.read_bytes() + b"\0")
    assert recorded_sha256(csv) != table_sha256(csv)


def test_general_expression_is_the_same_from_parts_and_a_single_file(tmp_path):
    (tmp_path / "parts").mkdir()
    (tmp_path / "single").mkdir()
    parts_csv = write_csv(tmp_path / "parts", cli.GENE_EXPRESSION_CSV, GENES)
    single_csv = write_csv(tmp_path / "single", cli.GENE_EXPRESSION_CSV, GENES)
    write_table(parts_csv, "pandas")
    pd.read_csv(single_csv).to_parquet(single_csv.with_suffix(".parquet"), index=False)

    valid = {"TSPAN6", "DPM1", "SCYL3", "FGR"}
    got = get_general_expression_of_genes(parts_csv, valid)
    want = get_general_expression_of_genes(single_csv, valid)
    pd.testing.assert_frame_equal(got.reset_index(drop=True), want.reset_index(drop=True))
    assert set(got.Gene) == {g for g in GENES if g.split(" (")[0] in valid}


@pytest.fixture
def depmap_dir(tmp_path, monkeypatch):
    monkeypatch.setattr(cli, "get_data_dir", lambda: str(tmp_path))
    monkeypatch.setattr(cli, "_ensure_depmap_file", lambda *args: False)
    return tmp_path


def test_setup_depmap_writes_the_gene_table_once(depmap_dir):
    csv = write_csv(depmap_dir, cli.GENE_EXPRESSION_CSV, GENES)

    cli.setup_depmap.callback(force=False)
    parts = sorted((depmap_dir / (csv.stem + "_parts")).glob("part-*.parquet"))
    assert len(parts) == 4 and not csv.exists()
    assert recorded_sha256(csv) == table_sha256(csv)

    # A second run finds the table matching its hash, and neither downloads nor converts.
    before = {p: p.stat().st_mtime_ns for p in parts}
    cli.setup_depmap.callback(force=False)
    assert {p: p.stat().st_mtime_ns for p in parts} == before


def test_setup_depmap_keeps_the_single_file_of_an_older_data_dir(depmap_dir):
    csv = write_csv(depmap_dir, cli.GENE_EXPRESSION_CSV, GENES)
    single = csv.with_suffix(".parquet")
    pd.read_csv(csv).to_parquet(single, index=False)
    csv.unlink()
    # The hash as setup-depmap recorded it before the split: of the file itself.
    Path(f"{single}.sha256").write_text(sha256_file(str(single)))

    cli.setup_depmap.callback(force=False)

    assert single.exists() and not (depmap_dir / (csv.stem + "_parts")).exists()
    assert recorded_sha256(csv) == table_sha256(csv)
