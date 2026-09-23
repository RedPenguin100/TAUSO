"""A DepMap CSV converted by depmap_parquet reads back as the CSV, however the Parquet is stored."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from tauso import cli
from tauso.cli_utils import sha256_file
from tauso.expression import depmap_parquet
from tauso.expression.depmap_parquet import (
    csv_to_parquet,
    expression_columns,
    mean_expression,
    parquet_exists,
    parquet_sha256,
    read_cell_lines,
    read_columns,
    save_parquet_sha256,
    saved_parquet_sha256,
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
TABLES = {"transcript": (cli.TRANSCRIPT_EXPRESSION_CSV, TRANSCRIPTS), "gene": (cli.GENE_EXPRESSION_CSV, GENES)}


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


@pytest.fixture(autouse=True)
def small_parts(monkeypatch):
    monkeypatch.setattr(depmap_parquet, "COLUMNS_PER_PART", 2)


@pytest.mark.parametrize("table", TABLES)
def test_the_table_reads_back_as_the_csv(tmp_path, table):
    name, columns = TABLES[table]
    csv = write_csv(tmp_path, name, columns)
    whole = pd.read_csv(csv)

    csv_to_parquet(csv)
    ids, got_columns, values = read_cell_lines(csv, {"ACH-000001", "ACH-000002"})

    assert not csv.exists() and parquet_exists(csv)
    assert expression_columns(csv) == got_columns == columns
    assert list(ids) == ["ACH-000001", "ACH-000002"]
    # ACH-000002 is read from its default profile, the third row; the missing value reads as 0.
    want = whole.iloc[[0, 2]][columns].fillna(0.0).to_numpy(dtype=np.float64)
    np.testing.assert_array_equal(values, want)


@pytest.mark.parametrize("table", TABLES)
def test_column_means_are_in_the_order_asked(tmp_path, table):
    name, columns = TABLES[table]
    csv = write_csv(tmp_path, name, columns)
    whole = pd.read_csv(csv)
    csv_to_parquet(csv)

    asked = [columns[4], columns[0], columns[2]]
    np.testing.assert_array_equal(mean_expression(csv, asked), [np.nanmean(whole[c].to_numpy(float)) for c in asked])


@pytest.mark.parametrize("table", TABLES)
def test_read_columns_returns_every_row_as_asked(tmp_path, table):
    name, columns = TABLES[table]
    csv = write_csv(tmp_path, name, columns)
    whole = pd.read_csv(csv)
    csv_to_parquet(csv)

    asked = [columns[5], columns[2]]
    got = read_columns(csv, asked)

    # Every profile is a row, the non-default one included, and a missing value stays missing.
    assert list(got.ModelID) == ["ACH-000001", "ACH-000002", "ACH-000002", "ACH-000003"]
    assert list(got.columns[-2:]) == asked
    pd.testing.assert_frame_equal(got[asked], whole[asked])
    assert got[columns[2]].isna().sum() == 1
    with pytest.raises(KeyError):
        read_columns(csv, ["NOT A COLUMN"])


def test_an_older_single_file_reads_the_same(tmp_path):
    (tmp_path / "parts").mkdir()
    (tmp_path / "single").mkdir()
    parts_csv = write_csv(tmp_path / "parts", cli.GENE_EXPRESSION_CSV, GENES)
    single_csv = write_csv(tmp_path / "single", cli.GENE_EXPRESSION_CSV, GENES)
    csv_to_parquet(parts_csv)
    pd.read_csv(single_csv).to_parquet(single_csv.with_suffix(".parquet"), index=False)

    ids = {"ACH-000001", "ACH-000002", "ACH-000003"}
    got, want = read_cell_lines(parts_csv, ids), read_cell_lines(single_csv, ids)
    assert list(got[0]) == list(want[0]) and got[1] == want[1]
    np.testing.assert_array_equal(got[2], want[2])
    np.testing.assert_array_equal(mean_expression(parts_csv, GENES), mean_expression(single_csv, GENES))
    pd.testing.assert_frame_equal(read_columns(parts_csv, GENES[:3]), read_columns(single_csv, GENES[:3]))


def test_writing_replaces_an_older_single_file(tmp_path):
    csv = write_csv(tmp_path, cli.GENE_EXPRESSION_CSV, GENES)
    single = csv.with_suffix(".parquet")
    pd.read_csv(csv).to_parquet(single, index=False)
    save_parquet_sha256(csv)

    csv_to_parquet(csv)

    assert not single.exists() and not Path(f"{single}.sha256").exists()
    assert saved_parquet_sha256(csv) == parquet_sha256(csv)


def test_a_killed_conversion_is_not_reused(tmp_path):
    csv = write_csv(tmp_path, cli.TRANSCRIPT_EXPRESSION_CSV, TRANSCRIPTS)
    partial = tmp_path / (csv.stem + "_parts.partial")
    partial.mkdir()
    (partial / "part-000.parquet").write_text("half")
    assert not parquet_exists(csv)

    csv_to_parquet(csv)

    assert not partial.exists()
    assert read_cell_lines(csv, {"ACH-000003"})[1] == TRANSCRIPTS


def test_the_recorded_hash_notices_a_changed_part(tmp_path):
    csv = write_csv(tmp_path, cli.GENE_EXPRESSION_CSV, GENES)
    csv_to_parquet(csv)
    assert saved_parquet_sha256(csv) == parquet_sha256(csv)

    part = sorted((tmp_path / (csv.stem + "_parts")).glob("part-*.parquet"))[1]
    part.write_bytes(part.read_bytes() + b"\0")
    assert saved_parquet_sha256(csv) != parquet_sha256(csv)


def test_general_expression_is_the_same_from_parts_and_a_single_file(tmp_path):
    (tmp_path / "parts").mkdir()
    (tmp_path / "single").mkdir()
    parts_csv = write_csv(tmp_path / "parts", cli.GENE_EXPRESSION_CSV, GENES)
    single_csv = write_csv(tmp_path / "single", cli.GENE_EXPRESSION_CSV, GENES)
    csv_to_parquet(parts_csv)
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
    assert saved_parquet_sha256(csv) == parquet_sha256(csv)

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
    assert saved_parquet_sha256(csv) == parquet_sha256(csv)
