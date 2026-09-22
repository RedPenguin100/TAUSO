"""The transcript table split into column parts reads back as the same table."""

import numpy as np
import pyarrow.csv as pacsv
import pyarrow.parquet as pq
import pytest

from tauso import cli

HEADER = ",SequencingID,ModelID,IsDefaultEntryForModel,ModelConditionID,IsDefaultEntryForMC"
TRANSCRIPTS = [f"ENST{n:011d}.1" for n in range(7)]
ROWS = [
    ("0", "S1", "ACH-000001", "Yes", "MC1", "Yes"),
    ("1", "S2", "ACH-000002", "No", "MC2", "No"),
    ("2", "S3", "ACH-000002", "Yes", "MC3", "Yes"),
    ("3", "S4", "ACH-000003", "Yes", "MC4", "Yes"),
]


@pytest.fixture
def csv_path(tmp_path):
    rng = np.random.default_rng(0)
    lines = [HEADER + "," + ",".join(TRANSCRIPTS)]
    for profile in ROWS:
        values = [repr(float(v)) for v in rng.random(len(TRANSCRIPTS))]
        values[3] = ""  # a missing value reads as 0
        lines.append(",".join(profile) + "," + ",".join(values))
    path = tmp_path / cli.TRANSCRIPT_EXPRESSION_CSV
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture
def small_parts(monkeypatch):
    monkeypatch.setattr(cli, "TRANSCRIPT_COLUMNS_PER_PART", 3)
    monkeypatch.setattr(cli, "_ensure_depmap_file", lambda *args: False)


def test_parts_read_back_as_the_whole_table(tmp_path, csv_path, small_parts):
    whole = tmp_path / "whole.parquet"
    pq.write_table(pacsv.read_csv(csv_path), whole)
    want = cli._read_transcript_cohort([str(whole)], {"ACH-000001", "ACH-000002"})

    parts = cli._ensure_transcript_parquet(str(tmp_path))
    got = cli._read_transcript_cohort(parts, {"ACH-000001", "ACH-000002"})

    # One part of profile columns, then 7 transcripts in parts of 3.
    assert len(parts) == 4
    assert list(got[0]) == list(want[0]) == ["ACH-000001", "ACH-000002"]
    assert got[1] == want[1] == TRANSCRIPTS
    np.testing.assert_array_equal(got[2], want[2])
    assert (got[2][:, 3] == 0).all()


def test_the_default_profile_row_is_the_one_read(tmp_path, csv_path, small_parts):
    table = pacsv.read_csv(csv_path).to_pandas()
    parts = cli._ensure_transcript_parquet(str(tmp_path))
    ids, _, values = cli._read_transcript_cohort(parts, {"ACH-000002"})

    default = table[(table.ModelID == "ACH-000002") & (table.IsDefaultEntryForModel == "Yes")]
    assert list(ids) == ["ACH-000002"]
    np.testing.assert_array_equal(values[0], default[TRANSCRIPTS].fillna(0).to_numpy()[0])


def test_the_csv_is_dropped_and_a_killed_conversion_is_not_reused(tmp_path, csv_path, small_parts):
    (tmp_path / (cli.TRANSCRIPT_EXPRESSION_PARTS + ".partial")).mkdir()
    (tmp_path / (cli.TRANSCRIPT_EXPRESSION_PARTS + ".partial") / "part-000.parquet").write_text("half")

    parts = cli._ensure_transcript_parquet(str(tmp_path))

    assert not csv_path.exists()
    assert not (tmp_path / (cli.TRANSCRIPT_EXPRESSION_PARTS + ".partial")).exists()
    assert all(pq.read_schema(p) for p in parts)
