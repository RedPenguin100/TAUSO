"""
Unit coverage for sum_exp_energy_by_trigger_pyarrow — the pyarrow RIsearch-output
parser the single-target-gene path reaches via off_target_specific_seq_pandarallel
/ on_target_total_hybridization (stream=True) -> _score_one_gene_streaming.

Unlike min_energy_by_pair_pyarrow (which returns per-(trigger,target) MIN energy),
this returns {trigger: sum(exp(-rt*energy))} over ALL hits — the score
_sum_exp_energy_by_trigger produced from a pandas DataFrame.

Like the min tests, these mock subprocess.Popen to stream synthetic RIsearch TSV
(no binary). The sum is floating-point order-dependent, so equality with the
pandas oracle is asserted within a relative tolerance (not bit-for-bit).
"""

import io
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

import tauso.features.hybridization.fast_hybridization as fh
from tauso.features.hybridization.fast_hybridization import sum_exp_energy_by_trigger_pyarrow

_RT = 0.616  # mirrors off_target_specific_gene._RT

# (trigger, target, energy). The single-gene path groups by trigger only and
# sums exp over every hit; the target column is irrelevant to the aggregation.
# Multiple hits per trigger, spread out so a small block_size splits a trigger's
# hits across pyarrow blocks (exercises per-block-sum -> concat -> final-sum).
HITS = [
    ("0", "GENE_A", -10.5),
    ("1", "GENE_A", -6.0),
    ("0", "GENE_A", -3.0),
    ("2", "GENE_A", -1.0),
    ("1", "GENE_A", -3.2),
    ("0", "GENE_A", -8.0),
    ("1", "GENE_A", -9.9123),
    ("2", "GENE_A", -2.5),
    ("0", "GENE_A", -0.75),
    ("1", "GENE_A", -5.25),
]


def _tsv_bytes(hits) -> bytes:
    lines = [f"{t}\t1\t20\t{g}\t100\t120\t900\t{e}" for t, g, e in hits]
    return ("\n".join(lines) + "\n").encode() if lines else b""


def _patched_popen(stdout_bytes: bytes, returncode: int = 0):
    proc = MagicMock()
    proc.stdout = io.BytesIO(stdout_bytes)
    proc.wait.return_value = returncode
    cm = MagicMock()
    cm.__enter__.return_value = proc
    cm.__exit__.return_value = False
    return patch.object(fh.subprocess, "Popen", MagicMock(return_value=cm))


def _pandas_oracle(hits, rt=_RT) -> dict:
    """The semantics the pyarrow path must reproduce: numpy exp + pandas groupby sum."""
    df = pd.DataFrame(hits, columns=["trigger", "target", "energy"])
    df["trigger"] = df["trigger"].astype(str)
    exp_vals = np.exp(-rt * df["energy"].to_numpy())
    s = pd.Series(exp_vals).groupby(df["trigger"].to_numpy(), sort=False).sum()
    return {str(k): float(v) for k, v in s.items()}


def _run(hits, block_size):
    with _patched_popen(_tsv_bytes(hits)):
        return sum_exp_energy_by_trigger_pyarrow(
            trigger_id_seq_pairs=[("0", "ACGTACGTACGTACGTACGT")],
            target_file_path="FAKE",
            block_size=block_size,
            transpose=True,
            rt=_RT,
        )


def _assert_close(result, expected, rel=1e-9):
    assert set(result) == set(expected)
    for k in expected:
        assert result[k] == pytest.approx(expected[k], rel=rel), k


# block_size 128 forces several blocks over these rows; 64<<20 is a single block.
@pytest.mark.parametrize("block_size", [128, 64 << 20])
def test_pyarrow_matches_pandas_sum_exp(block_size):
    _assert_close(_run(HITS, block_size), _pandas_oracle(HITS))


def test_multiblock_equals_single_block():
    """A trigger's hits split across blocks must sum to the same as one big block."""
    multi = _run(HITS, block_size=128)
    single = _run(HITS, block_size=64 << 20)
    _assert_close(multi, single, rel=1e-12)


def test_empty_output_returns_empty_dict():
    with _patched_popen(b""):
        result = sum_exp_energy_by_trigger_pyarrow(
            trigger_id_seq_pairs=[("0", "ACGTACGTACGTACGTACGT")],
            target_file_path="FAKE",
        )
    assert result == {}


def test_empty_pairs_short_circuits():
    assert sum_exp_energy_by_trigger_pyarrow(trigger_id_seq_pairs=[], target_file_path="FAKE") == {}
