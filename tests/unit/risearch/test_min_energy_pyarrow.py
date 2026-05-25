"""
Unit coverage for min_energy_by_pair_pyarrow — the pyarrow RIsearch-output parser
that production reaches via populate_off_target_* (stream=True).

The existing risearch unit tests exercise the non-streaming compute_group_batch
(stream=False), which production never runs, so the pyarrow path had no unit
coverage. These tests close that gap without invoking the RIsearch binary by
mocking subprocess.Popen to stream synthetic TSV from a BytesIO.

What they pin down:
  - the {(trigger, target): min_energy} result is bit-for-bit equal to a pandas
    groupby(...).min() on the same bytes (the "float parse matches pandas, MIN is
    exact" claim in the docstring);
  - the multi-block concat+regroup branch (the whole reason pyarrow exists for
    cutoff=0's millions of rows) is correct — i.e. a pair whose strongest hit
    lands in a later pyarrow block than an earlier weaker hit still gets the
    global min;
  - empty RIsearch stdout returns {}.
"""

import io
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

import tauso.features.hybridization.fast_hybridization as fh
from tauso.features.hybridization.fast_hybridization import min_energy_by_pair_pyarrow

COLUMNS = [
    "trigger",
    "trigger_start",
    "trigger_end",
    "target",
    "target_start",
    "target_end",
    "score",
    "energy",
]

# (trigger, target, energy) hits. Multiple hits per pair, and the *strongest*
# hit for a pair is deliberately placed far from its weaker hits so that with a
# small block_size they fall in different pyarrow blocks — exercising the
# per-block-min -> concat -> final-min path.
HITS = [
    ("0", "GENE_A", -3.0),
    ("0", "GENE_B", -5.0),
    ("1", "GENE_C", -3.2),
    ("0", "GENE_A", -10.5),  # min for (0, GENE_A), placed later
    ("2", "GENE_B", -1.0),
    ("1", "GENE_A", -6.0),
    ("0", "GENE_B", -7.5),  # min for (0, GENE_B), placed later
    ("1", "GENE_C", -9.9123),  # min for (1, GENE_C), several decimals
    ("2", "GENE_B", -2.5),
    ("0", "GENE_A", -8.0),
]


def _tsv_bytes(hits) -> bytes:
    lines = []
    for trigger, target, energy in hits:
        lines.append(f"{trigger}\t1\t20\t{target}\t100\t120\t900\t{energy}")
    return ("\n".join(lines) + "\n").encode() if lines else b""


def _patched_popen(stdout_bytes: bytes, returncode: int = 0):
    """A context-manager mock matching how min_energy_by_pair_pyarrow uses Popen."""
    proc = MagicMock()
    proc.stdout = io.BytesIO(stdout_bytes)
    proc.wait.return_value = returncode
    cm = MagicMock()
    cm.__enter__.return_value = proc
    cm.__exit__.return_value = False
    popen = MagicMock(return_value=cm)
    return patch.object(fh.subprocess, "Popen", popen)


def _pandas_oracle(hits) -> dict:
    """The semantics the pyarrow path must reproduce: pandas groupby min."""
    df = pd.DataFrame(hits, columns=["trigger", "target", "energy"])
    df["trigger"] = df["trigger"].astype(str)
    df["target"] = df["target"].astype(str)
    s = df.groupby(["trigger", "target"], sort=False)["energy"].min()
    return {k: float(v) for k, v in s.items()}


def _run(hits, block_size):
    with _patched_popen(_tsv_bytes(hits)):
        return min_energy_by_pair_pyarrow(
            trigger_id_seq_pairs=[("0", "ACGTACGTACGTACGTACGT")],
            target_file_path="FAKE",
            block_size=block_size,
            transpose=True,
        )


# block_size 128 forces ~13 blocks over these rows; 64<<20 is a single block.
@pytest.mark.parametrize("block_size", [128, 64 << 20])
def test_pyarrow_matches_pandas_groupby_min(block_size):
    result = _run(HITS, block_size)
    expected = _pandas_oracle(HITS)
    assert result == expected  # exact: float parse + MIN must match pandas bit-for-bit


def test_multiblock_equals_single_block():
    """Cross-block min: tiny blocks (concat path) must equal one big block."""
    multi = _run(HITS, block_size=128)
    single = _run(HITS, block_size=64 << 20)
    assert multi == single
    # sanity: the deliberately-late strongest hits won as the min
    assert multi[("0", "GENE_A")] == -10.5
    assert multi[("0", "GENE_B")] == -7.5
    assert multi[("1", "GENE_C")] == -9.9123


def test_empty_output_returns_empty_dict():
    with _patched_popen(b""):
        result = min_energy_by_pair_pyarrow(
            trigger_id_seq_pairs=[("0", "ACGTACGTACGTACGTACGT")],
            target_file_path="FAKE",
        )
    assert result == {}


def test_empty_pairs_short_circuits():
    assert min_energy_by_pair_pyarrow(trigger_id_seq_pairs=[], target_file_path="FAKE") == {}
