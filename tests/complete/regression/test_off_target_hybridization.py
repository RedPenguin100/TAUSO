import pytest

from tauso.features.hybridization.off_target.off_target_specific_gene import (
    off_target_single_gene_hybridization,
    on_target_total_hybridization,
)
from tauso.features.hybridization.off_target.on_target_multiplicity import (
    on_target_log_number_of_sites,
)
from tests.complete.conftest import get_n_jobs

SINGLE_TARGET_GENES = ["RNASEH1"]
CUTOFFS = [800, 1000, 1200]
# Production off-target cutoffs; see calculator.calculate_off_target_single.
SINGLE_OFF_TARGET_CUTOFFS = [800, 1000, 1200]
RRNA_CUTOFFS = [800, 1000, 1200]


@pytest.mark.parametrize("mini_structure_data", [1000], indirect=True)
def test_on_target_hybridization_regression(mini_structure_data, gene_to_data, dataframe_regression):
    data = mini_structure_data.copy()
    data, feature_names = on_target_total_hybridization(data, gene_to_data, cutoffs=CUTOFFS, n_jobs=get_n_jobs())
    # sum(exp(-energy/RT)) is float-order dependent under threads (~1e-5)
    dataframe_regression.check(data[["index_oligo"] + feature_names], default_tolerance={"atol": 1e-4, "rtol": 1e-4})


@pytest.mark.parametrize("mini_structure_data", [1000], indirect=True)
def test_on_target_multiplicity_regression(mini_structure_data, gene_to_data, dataframe_regression):
    """Site-resolved on-target multiplicity: log effective number of on-target sites per cutoff
    (the on_target_log_number_of_sites feature, backed by stats_by_trigger_multi_cutoff)."""
    data = mini_structure_data.copy()
    data, feature_names = on_target_log_number_of_sites(data, gene_to_data, cutoffs=CUTOFFS, n_jobs=get_n_jobs())
    # log_eff derives from sum(exp(-energy/RT)); float-order dependent under threads (~1e-5)
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
