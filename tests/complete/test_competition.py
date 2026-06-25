"""Regression tests for the competitor-score generators on a tiny sample.

Each test needs its external tool (OligoWalk/miranda/sfold binaries, or the 'pfred' Docker
container) and is skipped if it is unavailable -- so CI without the tools just skips them.
Regenerate baselines locally where the tools exist:  pytest tests/complete/test_competition.py --force-regen
"""
import pytest

from tauso.data.consts import ASO_SEQUENCE
from tests.complete.competitor_tools import pfred_container_running, tool_missing


@pytest.mark.skipif(tool_missing("OligoWalk"), reason="OligoWalk (rnastructure) not on PATH")
def test_oligowalk_regression(common_gene_sample, gene_to_data, dataframe_regression):
    from tauso.features.competition.oligowalk_utils import populate_oligowalk

    df, feats = populate_oligowalk(common_gene_sample.copy(), gene_to_data, num_workers=2)
    cols = ["index_oligo"] + [f for f in feats if f != "error"]
    dataframe_regression.check(df[cols].reset_index(drop=True))


@pytest.mark.skipif(tool_missing("miranda"), reason="miranda not on PATH")
def test_miranda_regression(common_gene_sample, gene_to_data, dataframe_regression):
    from notebooks.competitors.miranda_utils import get_populated_df_with_miranda_features

    df = get_populated_df_with_miranda_features(common_gene_sample.copy(), gene_to_data, max_workers=2)
    dataframe_regression.check(df[["index_oligo", "miranda_score", "miranda_energy"]].reset_index(drop=True))


@pytest.mark.skipif(tool_missing("sfold"), reason="sfold not on PATH")
def test_sfold_runs(common_gene_sample, gene_to_data):
    # sfold is stochastic (Boltzmann sampling), so sanity-check the output instead of pinning exact values.
    from notebooks.competitors.sfold_utils import calculate_sfold_accessibility

    acc = calculate_sfold_accessibility(df=common_gene_sample.copy(), gene_to_data=gene_to_data, n_cores=2)["sfold_accessibility"]
    assert len(acc) == len(common_gene_sample)
    assert acc.notna().any(), "sfold produced all-NaN -- the binary likely failed"
    valid = acc.dropna()
    assert ((valid >= 0) & (valid <= 1)).all(), "sfold_accessibility must be a probability in [0, 1]"


@pytest.mark.skipif(not pfred_container_running(), reason="'pfred' Docker container not running")
def test_pfred_regression(common_gene_sample, dataframe_regression, tmp_path):
    from notebooks.competitors.PFRED.pfred_glue import populate_pfred

    df, feats = populate_pfred(common_gene_sample.copy(), seq_col=ASO_SEQUENCE, temp_dir=str(tmp_path))
    dataframe_regression.check(df[["index_oligo"] + list(feats)].reset_index(drop=True))
