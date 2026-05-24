"""
Test compute_group_batch logic using synthetic RIsearch output (no binary invoked).
Covers:
  - O(N×M) fix: batch groupby aggregate + dict lookup
  - Self-target filtering per row
  - prebuilt_target_path skips file create/delete
  - ARTM score computed correctly
  - Empty group returns empty Series
"""

from unittest.mock import patch

import pandas as pd
import pytest

from tauso.data.consts import CANONICAL_GENE, SEQUENCE
from tauso.features.hybridization.off_target.add_off_target_feat import (
    AggregationMethod,
    compute_group_batch,
)

# Minimal fake RIsearch TSV output (trigger, t_start, t_end, target, ta_start, ta_end, score, energy)
FAKE_RISEARCH_OUTPUT = """\
0\t1\t20\tGENE_A\t100\t120\t1200\t-10.5
0\t1\t20\tGENE_B\t200\t220\t900\t-5.0
0\t1\t20\tGENE_A\t150\t170\t800\t-8.0
1\t1\t20\tGENE_B\t200\t220\t1100\t-7.5
1\t1\t20\tGENE_C\t50\t70\t950\t-3.2
2\t1\t20\tGENE_A\t100\t120\t1000\t-6.0
"""

# Expression map: {gene: (TPM, norm)}
EXP_MAP = {
    "GENE_A": (2_000_000, 0.5),  # TPM=2e6 → scaled = 2.0
    "GENE_B": (1_000_000, 0.3),  # TPM=1e6 → scaled = 1.0
    "GENE_C": (500_000, 0.2),  # TPM=5e5 → scaled = 0.5
}


def _make_group_df(rows):
    """rows: list of (canonical_gene, sequence)"""
    df = pd.DataFrame(rows, columns=[CANONICAL_GENE, SEQUENCE])
    return df


def _mock_batch(fake_output, prebuilt_target_path="FAKE"):
    """Patch get_triggers_mfe_scores_batch to return synthetic output, then call compute_group_batch."""
    with patch(
        "tauso.features.hybridization.off_target.add_off_target_feat.get_triggers_mfe_scores_batch",
        return_value=fake_output,
    ):
        return compute_group_batch


def test_artm_scores_computed_correctly():
    """ARTM: score = sum(energy * expr_tpm/1e6) for hits with energy<0 and expr>0, excluding own gene."""
    # Row 0: GENE_X (own). Hits to GENE_A (two hits → min energy -10.5), GENE_B (-5.0).
    # Row 1: GENE_A (own). Hits to GENE_B (-7.5), GENE_C (-3.2). Own gene GENE_A excluded.
    # Row 2: GENE_B (own). Hit to GENE_A (-6.0). Own gene GENE_B excluded.
    group_df = _make_group_df(
        [
            ("GENE_X", "AAAAAAAAAAAAAAAAAAAA"),
            ("GENE_A", "TTTTTTTTTTTTTTTTTTTT"),
            ("GENE_B", "GGGGGGGGGGGGGGGGGGGG"),
        ]
    )

    batch_fn = _mock_batch(FAKE_RISEARCH_OUTPUT)
    with patch(
        "tauso.features.hybridization.off_target.add_off_target_feat.get_triggers_mfe_scores_batch",
        return_value=FAKE_RISEARCH_OUTPUT,
    ):
        scores = compute_group_batch(
            group_df,
            seq_map={},
            exp_map=EXP_MAP,
            cutoff=0,
            method=AggregationMethod.ARTM,
            prebuilt_target_path="FAKE",
        )

    # Row 0 (GENE_X own): GENE_A min_energy=-10.5, GENE_B=-5.0. Neither is own gene.
    # ARTM: (-10.5 * 2.0) + (-5.0 * 1.0) = -21.0 + -5.0 = -26.0
    assert scores[0] == pytest.approx(-26.0, rel=1e-6), f"Row 0: {scores[0]}"

    # Row 1 (GENE_A own): GENE_B=-7.5, GENE_C=-3.2. GENE_A excluded.
    # ARTM: (-7.5 * 1.0) + (-3.2 * 0.5) = -7.5 + -1.6 = -9.1
    assert scores[1] == pytest.approx(-9.1, rel=1e-6), f"Row 1: {scores[1]}"

    # Row 2 (GENE_B own): GENE_A=-6.0. GENE_B excluded (no hit to GENE_B from row 2).
    # ARTM: (-6.0 * 2.0) = -12.0
    assert scores[2] == pytest.approx(-12.0, rel=1e-6), f"Row 2: {scores[2]}"


def test_self_target_always_excluded():
    """When every hit is to the row's own canonical gene, score must be 0."""
    group_df = _make_group_df(
        [
            ("GENE_A", "AAAAAAAAAAAAAAAAAAAA"),
        ]
    )
    # Only hit for row 0 is to GENE_A (its own gene)
    only_self_output = "0\t1\t20\tGENE_A\t100\t120\t1000\t-9.0\n"

    with patch(
        "tauso.features.hybridization.off_target.add_off_target_feat.get_triggers_mfe_scores_batch",
        return_value=only_self_output,
    ):
        scores = compute_group_batch(
            group_df,
            seq_map={},
            exp_map=EXP_MAP,
            cutoff=0,
            method=AggregationMethod.ARTM,
            prebuilt_target_path="FAKE",
        )

    assert scores[0] == pytest.approx(0.0)


def test_empty_risearch_output_returns_zeros():
    group_df = _make_group_df(
        [
            ("GENE_A", "AAAAAAAAAAAAAAAAAAAA"),
            ("GENE_B", "TTTTTTTTTTTTTTTTTTTT"),
        ]
    )
    with patch(
        "tauso.features.hybridization.off_target.add_off_target_feat.get_triggers_mfe_scores_batch",
        return_value="",
    ):
        scores = compute_group_batch(
            group_df,
            seq_map={},
            exp_map=EXP_MAP,
            cutoff=0,
            method=AggregationMethod.ARTM,
            prebuilt_target_path="FAKE",
        )

    assert (scores == 0.0).all()


def test_empty_group_df_returns_empty_series():
    group_df = pd.DataFrame(columns=[CANONICAL_GENE, SEQUENCE])
    scores = compute_group_batch(group_df, {}, {}, 0, AggregationMethod.ARTM)
    assert scores.empty


def test_prebuilt_target_path_skips_file_ops(tmp_path):
    """When prebuilt_target_path is provided, no target file is created or deleted."""
    group_df = _make_group_df([("GENE_X", "AAAAAAAAAAAAAAAAAAAA")])
    sentinel = tmp_path / "should_not_be_deleted.fa"
    sentinel.write_text(">dummy\nACGT\n")

    with patch(
        "tauso.features.hybridization.off_target.add_off_target_feat.get_triggers_mfe_scores_batch",
        return_value="",
    ), patch("tauso.features.hybridization.off_target.add_off_target_feat.dump_target_file") as mock_dump:
        compute_group_batch(
            group_df,
            seq_map={},
            exp_map={},
            cutoff=0,
            method=AggregationMethod.ARTM,
            prebuilt_target_path=str(sentinel),
        )
        mock_dump.assert_not_called()

    # The pre-built file must still exist (we don't own it)
    assert sentinel.exists()


def test_min_energy_per_target_used():
    """When multiple hits to the same target exist for one query, only the min energy is used."""
    group_df = _make_group_df([("GENE_X", "AAAAAAAAAAAAAAAAAAAA")])
    # Row 0 has two hits to GENE_A: -3.0 and -15.0 → min = -15.0
    multi_hit_output = "0\t1\t20\tGENE_A\t100\t120\t1000\t-3.0\n0\t1\t20\tGENE_A\t200\t220\t900\t-15.0\n"
    with patch(
        "tauso.features.hybridization.off_target.add_off_target_feat.get_triggers_mfe_scores_batch",
        return_value=multi_hit_output,
    ):
        scores = compute_group_batch(
            group_df,
            seq_map={},
            exp_map=EXP_MAP,
            cutoff=0,
            method=AggregationMethod.ARTM,
            prebuilt_target_path="FAKE",
        )

    # ARTM: -15.0 * (2_000_000 / 1e6) = -15.0 * 2.0 = -30.0
    assert scores[0] == pytest.approx(-30.0, rel=1e-6)
