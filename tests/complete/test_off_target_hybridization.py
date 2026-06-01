import pytest

from tauso.features.hybridization.off_target.add_off_target_feat import AggregationMethod
from tauso.features.hybridization.off_target.off_target_specific_gene import (
    off_target_single_gene_hybridization,
    on_target_total_hybridization,
)
from tauso.populate.populate_off_target import (
    populate_off_target_general,
    populate_off_target_specific,
)
from tests.complete.conftest import get_n_jobs

SINGLE_TARGET_GENES = ["RNASEH1"]
CUTOFFS = [800, 1000, 1200]
# Production off-target cutoffs; see calculator.calculate_off_target_single.
SINGLE_OFF_TARGET_CUTOFFS = [800, 1000, 1200]
RRNA_CUTOFFS = [800, 1000, 1200]


@pytest.fixture(scope="session")
def transcriptomes_with_general(transcriptomes):
    """Extends the transcriptomes dict with a 'general' key for use in general off-target tests."""
    from pathlib import Path

    from tauso.common.gtf import filter_gtf_genes
    from tauso.data.data import get_data_dir, load_gtf_db
    from tauso.features.expression.general_expression import get_general_expression_of_genes

    db = load_gtf_db()
    valid_genes = filter_gtf_genes(db, filter_mode="non_mt")

    data_dir = get_data_dir()
    exp_path = Path(data_dir) / "OmicsExpressionTPMLogp1HumanAllGenesStranded.csv"

    mean_exp_data = get_general_expression_of_genes(exp_path, valid_genes)

    result = dict(transcriptomes)
    result["general"] = mean_exp_data
    return result


@pytest.fixture(scope="session")
def mini_structure_data(request, structure_data):
    n_samples = getattr(request, "param", 1000)
    actual_samples = min(n_samples, len(structure_data))
    return structure_data.sample(n=actual_samples, random_state=42).copy()


@pytest.mark.parametrize("mini_structure_data", [1000], indirect=True)
def test_on_target_hybridization_regression(mini_structure_data, gene_to_data, dataframe_regression):
    data = mini_structure_data.copy()
    data, feature_names = on_target_total_hybridization(data, gene_to_data, cutoffs=CUTOFFS, n_jobs=get_n_jobs())
    # sum(exp(-RT·energy)) is float-order dependent under threads (~1e-5)
    dataframe_regression.check(data[["index_oligo"] + feature_names], default_tolerance={"atol": 1e-4, "rtol": 1e-4})


@pytest.mark.parametrize("mini_structure_data", [1000], indirect=True)
def test_off_target_single_regression(mini_structure_data, gene_to_data_full, dataframe_regression):
    data = mini_structure_data.copy()
    feature_names = []
    for target_gene in SINGLE_TARGET_GENES:
        data, fnames = off_target_single_gene_hybridization(
            data, target_gene, gene_to_data_full, cutoffs=SINGLE_OFF_TARGET_CUTOFFS, n_jobs=get_n_jobs()
        )
        feature_names += fnames
    dataframe_regression.check(data[["index_oligo"] + feature_names], default_tolerance={"atol": 1e-4, "rtol": 1e-4})


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


@pytest.mark.parametrize("mini_structure_data", [1000], indirect=True)
def test_off_target_single_rrna_regression(mini_structure_data, dataframe_regression):
    """rRNA off-target features + aggregated total across the shipped cutoffs. Skips if FASTA absent."""
    from tauso.features.hybridization.off_target.rrna_targets import RRNA_ACCESSIONS, get_rrna_loci

    try:
        rrna_loci = get_rrna_loci()
    except FileNotFoundError:
        pytest.skip("rRNA reference FASTA not present; run rrna_targets fetch to enable.")

    data = mini_structure_data.copy()
    rrna_species = list(RRNA_ACCESSIONS)
    feature_names = []
    for target_gene in rrna_species:
        data, fnames = off_target_single_gene_hybridization(
            data, target_gene, rrna_loci, cutoffs=RRNA_CUTOFFS, n_jobs=get_n_jobs()
        )
        feature_names += fnames
    for cutoff in RRNA_CUTOFFS:
        total_col = f"off_target_single_rRNA_total_c{cutoff}"
        data[total_col] = data[[f"off_target_single_{sp}_c{cutoff}" for sp in rrna_species]].sum(axis=1)
        feature_names.append(total_col)
    # FP-tolerant (parallel aggregation / parser wobble ~1e-5), as in test_off_target_specific_regression.
    dataframe_regression.check(data[["index_oligo"] + feature_names], default_tolerance={"atol": 1e-4, "rtol": 1e-4})


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
        method=AggregationMethod.ARTM,
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
        method=AggregationMethod.ARTM,
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
        method=AggregationMethod.ARTM,
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
        method=AggregationMethod.ARTM,
        n_jobs=get_n_jobs(),
    )
    multi, feats = populate_off_target_specific(ASO_df=mini_structure_data.copy(), cutoff_list=cutoffs, **base)
    for c in cutoffs:
        single, sfeats = populate_off_target_specific(ASO_df=mini_structure_data.copy(), cutoff_list=[c], **base)
        pd.testing.assert_series_equal(multi[sfeats[0]], single[sfeats[0]])
