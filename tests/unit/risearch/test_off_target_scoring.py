"""
Off-target scoring/aggregation unit tests (no RIsearch binary — synthetic inputs).

Targets the real production functions:
  - score_min_energy_per_cutoff: self-target filtering + ARTM scoring
  - min_energy_by_pair_multi_cutoff: min-energy per (trigger, target), cutoff filter on score
  - sum_exp_by_trigger_multi_cutoff: Boltzmann sum(exp(-RT*energy)) per trigger, cutoff filter

(Replaces the old test_compute_group_batch.py + test_scoring_vectorization.py, whose
target functions — compute_group_batch / _sum_exp_energy_by_trigger — were removed.)
"""

import numpy as np
import pyarrow as pa
import pytest

from tauso.features.hybridization.fast_hybridization import (
    min_energy_by_pair_multi_cutoff,
    sum_exp_by_trigger_multi_cutoff,
)
from tauso.features.hybridization.off_target.add_off_target_feat import (
    AggregationMethod,
    score_min_energy_per_cutoff,
)

_RT = 0.616

# (trigger, target, score, energy) — two hits to (0, GENE_A); a score==1000 boundary hit.
HITS = {
    "trigger": ["0", "0", "0", "1", "1", "2"],
    "target": ["GENE_A", "GENE_B", "GENE_A", "GENE_B", "GENE_C", "GENE_A"],
    "score": [1200, 900, 800, 1100, 950, 1000],
    "energy": [-10.5, -5.0, -8.0, -7.5, -3.2, -6.0],
}
EXP_MAP = {"GENE_A": (2_000_000, 0.5), "GENE_B": (1_000_000, 0.3), "GENE_C": (500_000, 0.2)}


def _batch():
    return pa.record_batch(
        {
            "trigger": pa.array(HITS["trigger"], pa.string()),
            "target": pa.array(HITS["target"], pa.string()),
            "score": pa.array(HITS["score"], pa.int64()),
            "energy": pa.array(HITS["energy"], pa.float64()),
        }
    )


def _reduce(agg):
    return agg.finalize(agg.combine(_batch()))


# --- min_energy_by_pair_multi_cutoff -----------------------------------------
def test_min_energy_keeps_strongest_per_pair():
    per_cutoff = _reduce(min_energy_by_pair_multi_cutoff([0]))
    # (0, GENE_A) had -10.5 and -8.0 -> keeps the min (-10.5)
    assert per_cutoff[0][("0", "GENE_A")] == pytest.approx(-10.5)
    assert per_cutoff[0][("0", "GENE_B")] == pytest.approx(-5.0)


def test_cutoff_filter_is_on_score_and_exclusive():
    per_cutoff = _reduce(min_energy_by_pair_multi_cutoff([0, 1000]))
    # score > 1000 keeps only (0,GENE_A,1200) and (1,GENE_B,1100); the score==1000 hit is excluded
    assert set(per_cutoff[1000]) == {("0", "GENE_A"), ("1", "GENE_B")}
    assert ("2", "GENE_A") in per_cutoff[0] and ("2", "GENE_A") not in per_cutoff[1000]


# --- sum_exp_by_trigger_multi_cutoff -----------------------------------------
def test_sum_exp_matches_reference():
    per_cutoff = _reduce(sum_exp_by_trigger_multi_cutoff([0], rt=_RT))
    energies = np.array(HITS["energy"])
    triggers = np.array(HITS["trigger"])
    for t in ("0", "1", "2"):
        expected = float(np.exp(-_RT * energies[triggers == t]).sum())
        assert per_cutoff[0][t] == pytest.approx(expected, rel=1e-9)


# --- score_min_energy_per_cutoff (self-filter + ARTM) ------------------------
def test_artm_score_excludes_own_gene():
    """ARTM: sum(energy * expr_TPM/1e6) over hits with energy<0, expr>0, excluding own gene."""
    per_cutoff = {
        0: {
            ("0", "GENE_A"): -10.5,
            ("0", "GENE_B"): -5.0,
            ("1", "GENE_B"): -7.5,
            ("1", "GENE_C"): -3.2,
            ("2", "GENE_A"): -6.0,
        }
    }
    indices = [0, 1, 2]
    idx_to_gene = {"0": "GENE_X", "1": "GENE_A", "2": "GENE_B"}  # own genes
    out = score_min_energy_per_cutoff(per_cutoff, indices, idx_to_gene, EXP_MAP, AggregationMethod.ARTM, [0])
    s = out[0]
    assert s[0] == pytest.approx(-10.5 * 2.0 + -5.0 * 1.0)  # GENE_X owns nothing here
    assert s[1] == pytest.approx(-7.5 * 1.0 + -3.2 * 0.5)  # own GENE_A would be excluded
    assert s[2] == pytest.approx(-6.0 * 2.0)


def test_self_hit_excluded():
    per_cutoff = {0: {("0", "GENE_A"): -9.0, ("0", "GENE_B"): -4.0}}
    out = score_min_energy_per_cutoff(per_cutoff, [0], {"0": "GENE_A"}, EXP_MAP, AggregationMethod.ARTM, [0])
    assert out[0][0] == pytest.approx(-4.0 * 1.0)  # GENE_A is own -> only GENE_B contributes


def test_empty_returns_zeros():
    out = score_min_energy_per_cutoff({0: {}}, [0, 1], {"0": "G", "1": "G"}, EXP_MAP, AggregationMethod.ARTM, [0])
    assert out[0].tolist() == [0.0, 0.0]
