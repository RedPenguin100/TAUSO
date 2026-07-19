import numpy as np
import pytest

from tauso.features.hybridization.off_target.add_off_target_feat import AggregationMethod
from tauso.features.hybridization.off_target.off_target_specific_gene import (
    add_on_target_site_features,
    off_target_single_gene_hybridization,
    on_target_total_hybridization,
)
from tauso.features.hybridization.off_target.on_target_multiplicity import (
    on_target_log_number_of_sites,
)
from tauso.populate.populate_off_target import (
    populate_off_target_general,
    populate_off_target_specific,
)
from tests.complete.conftest import get_n_jobs

CUTOFFS = [800, 1000, 1200]


@pytest.mark.parametrize("mini_structure_data", [1000], indirect=True)
def test_add_on_target_site_features_matches_separate(mini_structure_data, gene_to_data):
    """The single-scan add_on_target_site_features emits the same total-hybridization and
    log-number-of-sites columns as running the two per-feature entry points separately."""
    n_jobs = get_n_jobs()
    both, both_names = add_on_target_site_features(
        mini_structure_data.copy(), gene_to_data, cutoffs=CUTOFFS, n_jobs=n_jobs
    )
    tot, tot_names = on_target_total_hybridization(
        mini_structure_data.copy(), gene_to_data, cutoffs=CUTOFFS, n_jobs=n_jobs
    )
    log, log_names = on_target_log_number_of_sites(
        mini_structure_data.copy(), gene_to_data, cutoffs=CUTOFFS, n_jobs=n_jobs
    )
    for name in tot_names:
        assert name in both_names
        np.testing.assert_allclose(both[name].to_numpy(), tot[name].to_numpy(), rtol=1e-4, atol=1e-4)
    for name in log_names:
        assert name in both_names
        np.testing.assert_allclose(both[name].to_numpy(), log[name].to_numpy(), rtol=1e-4, atol=1e-4)


@pytest.mark.parametrize("mini_structure_data", [300], indirect=True)
def test_on_target_multi_cutoff_matches_per_cutoff(mini_structure_data, gene_to_data):
    """on_target deriving all cutoffs from one loose pass equals running each separately
    (within FP rounding of the order-dependent exp-sum)."""
    import pandas as pd

    cutoffs = [800, 1000, 1200]
    multi, _ = on_target_total_hybridization(
        mini_structure_data.copy(), gene_to_data, cutoffs=cutoffs, n_jobs=get_n_jobs()
    )
    for c in cutoffs:
        single, sfeats = on_target_total_hybridization(
            mini_structure_data.copy(), gene_to_data, cutoffs=[c], n_jobs=get_n_jobs()
        )
        pd.testing.assert_series_equal(multi[sfeats[0]], single[sfeats[0]])


def test_rrna_off_target_binder_vs_random():
    """Fast guard: an ASO complementary to 18S scores >0, a random one ~0. Skips if FASTA absent."""
    import pandas as pd

    from tauso.data.consts import ASO_SEQUENCE
    from tauso.features.hybridization.off_target.rrna_targets import get_rrna_loci
    from tauso.util import get_antisense

    try:
        loci = get_rrna_loci()
    except FileNotFoundError:
        pytest.skip("rRNA reference FASTA not present; run rrna_targets fetch to enable.")

    seq_18s = loci["rRNA_18S"].full_mrna
    binder = get_antisense(seq_18s[200:220])  # antisense of an 18S 20-mer -> binds 18S
    df = pd.DataFrame({ASO_SEQUENCE: [binder, "ACGTACGTACGTACGTACGT"]})

    df, feats = off_target_single_gene_hybridization(df, "rRNA_18S", loci, cutoffs=[1200], n_jobs=1)
    feat = feats[0]
    assert df[feat].iloc[0] > 0.0, "designed 18S binder should hybridize to 18S"
    assert df[feat].iloc[1] == 0.0, "random oligo should not strongly hybridize to 18S"


@pytest.mark.parametrize("mini_structure_data", [1000], indirect=True)
def test_off_target_specific_regression(mini_structure_data, gene_to_data, transcriptomes, dataframe_regression):
    data = mini_structure_data.copy()
    data, feature_names = populate_off_target_specific(
        ASO_df=data,
        gene_to_data=gene_to_data,
        cell_line2data=transcriptomes,
        top_n_list=[50],
        cutoff_list=[800],
        method=AggregationMethod.BOLTZMANN_SUM,
        n_jobs=get_n_jobs(),
    )
    # Parallel thread aggregation can produce floating-point differences ~1e-5
    dataframe_regression.check(data[["index_oligo"] + feature_names], default_tolerance={"atol": 1e-4, "rtol": 1e-4})


@pytest.mark.parametrize("mini_structure_data", [1000], indirect=True)
def test_off_target_general_regression(
    mini_structure_data, gene_to_data_full, transcriptomes_with_general, dataframe_regression
):
    data = mini_structure_data.copy()
    data, feature_names = populate_off_target_general(
        ASO_df=data,
        gene_to_data=gene_to_data_full,
        cell_line2data=transcriptomes_with_general,
        top_n_list=[25],
        cutoff_list=[800],
        method=AggregationMethod.BOLTZMANN_SUM,
        n_jobs=get_n_jobs(),
    )
    dataframe_regression.check(data[["index_oligo"] + feature_names])


@pytest.mark.parametrize("mini_structure_data", [300], indirect=True)
def test_off_target_general_multi_cutoff_matches_per_cutoff(
    mini_structure_data, gene_to_data_full, transcriptomes_with_general
):
    """off_target_general multi-cutoff collapse test: deriving every cutoff from one loose
    RIsearch pass equals running each cutoff on its own."""
    import pandas as pd

    cutoffs = [800, 1000, 1200]
    kwargs = dict(
        gene_to_data=gene_to_data_full,
        cell_line2data=transcriptomes_with_general,
        top_n_list=[25],
        method=AggregationMethod.BOLTZMANN_SUM,
        n_jobs=get_n_jobs(),
    )

    multi, _ = populate_off_target_general(ASO_df=mini_structure_data.copy(), cutoff_list=cutoffs, **kwargs)
    for c in cutoffs:
        single, feats = populate_off_target_general(ASO_df=mini_structure_data.copy(), cutoff_list=[c], **kwargs)
        pd.testing.assert_series_equal(multi[feats[0]], single[feats[0]])


@pytest.mark.parametrize("mini_structure_data", [300], indirect=True)
def test_off_target_specific_multi_cutoff_matches_per_cutoff(mini_structure_data, gene_to_data, transcriptomes):
    """Specific path: deriving all cutoffs from one loose pass equals running each separately."""
    import pandas as pd

    cutoffs = [800, 1000, 1200]
    base = dict(
        gene_to_data=gene_to_data,
        cell_line2data=transcriptomes,
        top_n_list=[50],
        method=AggregationMethod.BOLTZMANN_SUM,
        n_jobs=get_n_jobs(),
    )
    multi, feats = populate_off_target_specific(ASO_df=mini_structure_data.copy(), cutoff_list=cutoffs, **base)
    for c in cutoffs:
        single, sfeats = populate_off_target_specific(ASO_df=mini_structure_data.copy(), cutoff_list=[c], **base)
        pd.testing.assert_series_equal(multi[sfeats[0]], single[sfeats[0]])


@pytest.mark.parametrize("mini_structure_data", [300], indirect=True)
def test_off_target_general_multi_topn_matches_per_topn(
    mini_structure_data, gene_to_data_full, transcriptomes_with_general
):
    """general path: deriving multiple top_n from a single max(top_n) RIsearch pass
    (collapsed) equals running each top_n separately (uncollapsed)."""
    import pandas as pd

    top_ns = [25, 50, 100]
    kwargs = dict(
        gene_to_data=gene_to_data_full,
        cell_line2data=transcriptomes_with_general,
        cutoff_list=[800],
        method=AggregationMethod.BOLTZMANN_SUM,
        n_jobs=get_n_jobs(),
    )
    multi, _ = populate_off_target_general(ASO_df=mini_structure_data.copy(), top_n_list=top_ns, **kwargs)
    for n in top_ns:
        single, feats = populate_off_target_general(ASO_df=mini_structure_data.copy(), top_n_list=[n], **kwargs)
        pd.testing.assert_series_equal(multi[feats[0]], single[feats[0]])


@pytest.mark.parametrize("mini_structure_data", [300], indirect=True)
def test_off_target_specific_multi_topn_matches_per_topn(mini_structure_data, gene_to_data, transcriptomes):
    """specific path: deriving multiple top_n from a single max(top_n) RIsearch pass
    (collapsed) equals running each top_n separately (uncollapsed)."""
    import pandas as pd

    top_ns = [25, 50, 100]
    base = dict(
        gene_to_data=gene_to_data,
        cell_line2data=transcriptomes,
        cutoff_list=[800],
        method=AggregationMethod.BOLTZMANN_SUM,
        n_jobs=get_n_jobs(),
    )
    multi, _ = populate_off_target_specific(ASO_df=mini_structure_data.copy(), top_n_list=top_ns, **base)
    for n in top_ns:
        single, feats = populate_off_target_specific(ASO_df=mini_structure_data.copy(), top_n_list=[n], **base)
        pd.testing.assert_series_equal(multi[feats[0]], single[feats[0]])
