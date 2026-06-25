"""Regression tests for the competitor-score generators on a tiny sample.

Each test needs its external tool (OligoWalk/miranda/sfold binaries, or the 'pfred' Docker
container) and is skipped if it is unavailable -- so CI without the tools just skips them.
Regenerate baselines locally where the tools exist:  pytest tests/complete/test_competition.py --force-regen
"""
import shutil
import subprocess

import pytest

from tauso.data.consts import ASO_SEQUENCE, CANONICAL_GENE_NAME
from tauso.genome.read_human_genome import get_locus_to_data_dict


@pytest.fixture(scope="module")
def comp_sample(base_data):
    """Tiny deterministic slice: a few ASOs from the 2 most-represented genes (fast tools + gene dict)."""
    genes = base_data[CANONICAL_GENE_NAME].value_counts().index[:2].tolist()
    sub = base_data[base_data[CANONICAL_GENE_NAME].isin(genes)]
    return sub.groupby(CANONICAL_GENE_NAME, group_keys=False).head(6).reset_index(drop=True)


@pytest.fixture(scope="module")
def comp_gene_to_data(comp_sample):
    return get_locus_to_data_dict(include_introns=True, gene_subset=comp_sample[CANONICAL_GENE_NAME].unique().tolist())


def _missing(tool):
    return shutil.which(tool) is None


def _pfred_running():
    try:
        r = subprocess.run(["docker", "inspect", "-f", "{{.State.Running}}", "pfred"], capture_output=True, text=True)
        return r.stdout.strip() == "true"
    except Exception:
        return False


@pytest.mark.skipif(_missing("OligoWalk"), reason="OligoWalk (rnastructure) not on PATH")
def test_oligowalk_regression(comp_sample, comp_gene_to_data, dataframe_regression):
    from tauso.features.competition.oligowalk_utils import populate_oligowalk

    df, feats = populate_oligowalk(comp_sample.copy(), comp_gene_to_data, num_workers=2)
    cols = ["index_oligo"] + [f for f in feats if f != "error"]
    dataframe_regression.check(df[cols].reset_index(drop=True))


@pytest.mark.skipif(_missing("miranda"), reason="miranda not on PATH")
def test_miranda_regression(comp_sample, comp_gene_to_data, dataframe_regression):
    from notebooks.competitors.miranda_utils import get_populated_df_with_miranda_features

    df = get_populated_df_with_miranda_features(comp_sample.copy(), comp_gene_to_data, max_workers=2)
    dataframe_regression.check(df[["index_oligo", "miranda_score", "miranda_energy"]].reset_index(drop=True))


@pytest.mark.skipif(_missing("sfold"), reason="sfold not on PATH")
def test_sfold_runs(comp_sample, comp_gene_to_data):
    # sfold is stochastic (Boltzmann sampling), so sanity-check the output instead of pinning exact values.
    from notebooks.competitors.sfold_utils import calculate_sfold_accessibility

    acc = calculate_sfold_accessibility(df=comp_sample.copy(), gene_to_data=comp_gene_to_data, n_cores=2)["sfold_accessibility"]
    assert len(acc) == len(comp_sample)
    assert acc.notna().any(), "sfold produced all-NaN -- the binary likely failed"
    valid = acc.dropna()
    assert ((valid >= 0) & (valid <= 1)).all(), "sfold_accessibility must be a probability in [0, 1]"


@pytest.mark.skipif(not _pfred_running(), reason="'pfred' Docker container not running")
def test_pfred_regression(comp_sample, dataframe_regression):
    from notebooks.competitors.PFRED.pfred_glue import populate_pfred

    df, feats = populate_pfred(comp_sample.copy(), seq_col=ASO_SEQUENCE)
    dataframe_regression.check(df[["index_oligo"] + list(feats)].reset_index(drop=True))
