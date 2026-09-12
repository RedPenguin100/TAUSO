"""Search failures must never become successful zero-hit or intergenic results."""

import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest
from click.testing import CliRunner

from tauso import cli
from tauso.off_target import search

SEQUENCE = "ACGTACGTACGT"


@pytest.fixture(params=["single", "bulk", "counts"])
def run_search(request, monkeypatch, tmp_path):
    monkeypatch.setattr(search, "get_bowtie_index_base", lambda **kwargs: "test-index")
    fasta = tmp_path / "queries.fasta"
    fasta.write_text(f">{SEQUENCE}\n{SEQUENCE}\n")
    runners = {
        "single": lambda: search.run_bowtie_search(SEQUENCE, max_mismatches=2),
        "bulk": lambda: search.run_bowtie_search_bulk(str(fasta), max_mismatches=2),
        "counts": lambda: search.count_offtarget_matches_bulk([SEQUENCE], max_mismatches=2),
    }
    empty_results = {
        "single": ([], {"mismatches0": 0, "mismatches1": 0, "mismatches2": 0}),
        "bulk": [],
        "counts": {SEQUENCE: {0: 0, 1: 0, 2: 0}},
    }
    return runners[request.param], empty_results[request.param]


def test_alignment_failure_raises_with_command_and_diagnostics(run_search, monkeypatch):
    errors = []

    def fail(cmd, **kwargs):
        error = subprocess.CalledProcessError(2, cmd, stderr="Could not read index")
        errors.append(error)
        raise error

    monkeypatch.setattr(search.subprocess, "run", fail)
    run, _ = run_search
    with pytest.raises(RuntimeError, match="Could not read index") as caught:
        run()
    assert caught.value.__cause__ is errors[0]
    assert str(errors[0]) in str(caught.value)


@pytest.mark.parametrize("stdout", ["", f"{SEQUENCE}\t4\t*\t0\t0\t*\t*\t0\t0\t{SEQUENCE}\t*\n"])
def test_successful_search_without_alignments_keeps_empty_result(run_search, monkeypatch, stdout):
    def succeed(cmd, **kwargs):
        assert kwargs["check"] is True
        if "-f" in cmd:
            Path(cmd[-1]).write_text(stdout)
        return subprocess.CompletedProcess(cmd, 0, stdout=stdout, stderr="")

    monkeypatch.setattr(search.subprocess, "run", succeed)
    run, expected = run_search
    assert run() == expected


@pytest.fixture
def hit():
    return {"chrom": "chr1", "start": 100, "end": 112, "strand": "-", "sequence": SEQUENCE, "mismatches": 0}


@pytest.mark.parametrize("stage", ["load", "query"])
def test_annotation_failures_propagate(monkeypatch, hit, stage):
    error = OSError("Cannot read annotation")

    def fail(*args, **kwargs):
        raise error

    if stage == "load":
        monkeypatch.setattr(search, "load_gene_intervals", fail)
    else:
        tree = SimpleNamespace(all_overlaps_both=fail)
        monkeypatch.setattr(search, "_annotation_index", lambda genome: ({}, {("chr1", "+"): tree}))
    with pytest.raises(OSError) as caught:
        search.annotate_hits([hit], genome="test-annotation-failure")
    assert caught.value is error


def test_successful_annotation_without_overlaps_is_intergenic(monkeypatch, hit):
    monkeypatch.setattr(search, "_annotation_index", lambda genome: ({}, {}))
    result = search.annotate_hits([hit])
    assert len(result) == 1
    assert result.iloc[0]["region_type"] == "Intergenic"
    assert result.iloc[0]["gene_name"] is None


def test_cli_reports_failed_search_without_writing_results(monkeypatch, tmp_path):
    index_dir = tmp_path / "GRCh38_bowtie_index"
    index_dir.mkdir()
    (index_dir / "SUCCESS").touch()
    monkeypatch.setattr(cli, "get_paths", lambda genome: {"fasta": str(tmp_path / "GRCh38.fa")})
    monkeypatch.setattr(search, "get_bowtie_index_base", lambda **kwargs: "test-index")

    def fail(cmd, **kwargs):
        raise subprocess.CalledProcessError(2, cmd, stderr="Could not read index")

    monkeypatch.setattr(search.subprocess, "run", fail)
    output = tmp_path / "hits.csv"
    result = CliRunner().invoke(cli.main, ["run-off-target", SEQUENCE, "--output", str(output)])
    assert result.exit_code == 1
    assert "Search failed:" in result.output
    assert "Could not read index" in result.output
    assert "No off-targets found" not in result.output
    assert "Search completed" not in result.output
    assert not output.exists()
