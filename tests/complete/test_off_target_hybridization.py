import pytest

from tests.complete.conftest import get_n_jobs
from tauso.features.hybridization.off_target.add_off_target_feat import AggregationMethod
from tauso.populate.populate_off_target import (
    populate_off_target_general,
    populate_off_target_specific,
)
from tauso.features.hybridization.off_target.off_target_specific_gene import (
    off_target_specific_seq_pandarallel,
    on_target_total_hybridization,
)

SINGLE_TARGET_GENES = ["RNASEH1", "ACTB"]
CUTOFFS = [0, 1200]
RRNA_CUTOFFS = [0, 800, 1000, 1200, 1500]


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
    feature_names = []
    for cutoff in CUTOFFS:
        data, feature_name = on_target_total_hybridization(data, gene_to_data, cutoff=cutoff, n_jobs=get_n_jobs())
        feature_names.append(feature_name)
    dataframe_regression.check(data[["index_oligo"] + feature_names])


@pytest.mark.parametrize("mini_structure_data", [1000], indirect=True)
def test_off_target_single_regression(mini_structure_data, gene_to_data_full, dataframe_regression):
    data = mini_structure_data.copy()
    feature_names = []
    for target_gene in SINGLE_TARGET_GENES:
        for cutoff in CUTOFFS:
            data, feature_name = off_target_specific_seq_pandarallel(
                data, target_gene, gene_to_data_full, cutoff=cutoff, n_jobs=get_n_jobs()
            )
            feature_names.append(feature_name)
    dataframe_regression.check(data[["index_oligo"] + feature_names])


@pytest.fixture(scope="session")
def mini_seq_data():
    """Lightweight 1000-oligo sample (sequence + index only) for rRNA scoring.

    rRNA scoring needs only the ASO sequence, so this avoids the heavy structure/genome
    fixtures (and the CPU they spend) used by the other regression tests.
    """
    import pandas as pd

    from notebooks.consts import OLIGO_CSV_INDEXED
    from notebooks.data.OligoAI.utility import standardize_oligo_ai_data

    df = standardize_oligo_ai_data(pd.read_csv(OLIGO_CSV_INDEXED))
    return df.sample(n=min(1000, len(df)), random_state=42).copy()


def test_off_target_single_rrna_regression(mini_seq_data, dataframe_regression):
    """Regression for the cytoplasmic rRNA off-target features + aggregated total across the
    shipped cutoffs. Skips if the rRNA reference FASTA hasn't been fetched into the data dir.
    """
    from tauso.features.hybridization.off_target.rrna_targets import RRNA_ACCESSIONS, get_rrna_loci

    try:
        rrna_loci = get_rrna_loci()
    except FileNotFoundError:
        pytest.skip("rRNA reference FASTA not present; run rrna_targets fetch to enable.")

    data = mini_seq_data.copy()
    rrna_species = list(RRNA_ACCESSIONS)
    feature_names = []
    for target_gene in rrna_species:
        for cutoff in RRNA_CUTOFFS:
            data, feature_name = off_target_specific_seq_pandarallel(
                data, target_gene, rrna_loci, cutoff=cutoff, n_jobs=get_n_jobs()
            )
            feature_names.append(feature_name)
    for cutoff in RRNA_CUTOFFS:
        total_col = f"off_target_single_rRNA_total_c{cutoff}"
        data[total_col] = data[[f"off_target_single_{sp}_c{cutoff}" for sp in rrna_species]].sum(axis=1)
        feature_names.append(total_col)
    # Tolerance mirrors test_off_target_specific_regression: parallel-thread aggregation
    # (and parser internals) can produce floating-point differences ~1e-5.
    dataframe_regression.check(
        data[["index_oligo"] + feature_names], default_tolerance={"atol": 1e-4, "rtol": 1e-4}
    )


def test_rrna_off_target_binder_vs_random():
    """Fast deterministic guard for the rRNA off-target wiring: an ASO complementary to
    18S scores > 0, a random one ~0. No heavy fixtures, runs in well under a second.

    Skips if the rRNA reference FASTA hasn't been fetched into the data dir.
    """
    import pandas as pd

    from tauso.data.consts import SEQUENCE
    from tauso.features.hybridization.off_target.rrna_targets import get_rrna_loci
    from tauso.util import get_antisense

    try:
        loci = get_rrna_loci()
    except FileNotFoundError:
        pytest.skip("rRNA reference FASTA not present; run rrna_targets fetch to enable.")

    seq_18s = loci["rRNA_18S"].full_mrna
    binder = get_antisense(seq_18s[200:220])  # antisense of an 18S 20-mer -> binds 18S
    df = pd.DataFrame({SEQUENCE: [binder, "ACGTACGTACGTACGTACGT"]})

    df, feat = off_target_specific_seq_pandarallel(df, "rRNA_18S", loci, cutoff=1200, n_jobs=1)
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
