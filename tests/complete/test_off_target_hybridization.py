import pytest

from tests.complete.conftest import get_n_jobs
from tauso.features.hybridization_off_target.add_off_target_feat import AggregationMethod
from tauso.features.hybridization_off_target.off_target_feature import (
    populate_off_target_general,
    populate_off_target_specific,
)
from tauso.features.hybridization_off_target.off_target_specific_gene import (
    off_target_specific_seq_pandarallel,
    on_target_total_hybridization,
)

SINGLE_TARGET_GENES = ["RNASEH1", "ACTB"]
CUTOFFS = [0, 1200]


@pytest.fixture(scope="session")
def transcriptomes_with_general(transcriptomes):
    """Extends the transcriptomes dict with a 'general' key for use in general off-target tests."""
    from pathlib import Path

    from tauso.common.gtf import filter_gtf_genes
    from tauso.data.data import get_data_dir, load_gtf_db
    from tauso.features.hybridization_off_target.common import get_general_expression_of_genes

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
