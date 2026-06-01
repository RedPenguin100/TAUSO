import pytest

from tauso.data.consts import (
    STRUCT_SENSE_IN_CDS,
    STRUCT_SENSE_IN_CDS_NON_EXCLUSIVE,
    STRUCT_SENSE_IN_EXON,
    STRUCT_SENSE_IN_INTRON,
    STRUCTURE_SENSE_LENGTH,
    STRUCTURE_SENSE_START,
    STRUCTURE_SENSE_START_FROM_END,
    STRUCTURE_SENSE_TYPE,
    STRUCT_SENSE_IN_UTR,
)
from tauso.populate.populate_structure import get_populated_df_with_structure_features
from tauso.timer import Timer


@pytest.fixture
def sampled_base_data(request, base_data):
    """Samples the base data right before the test runs."""
    n_samples = getattr(request, "param", 1000)
    actual_samples = min(n_samples, len(base_data))
    return base_data.sample(n=actual_samples, random_state=42).copy()


@pytest.mark.parametrize("sampled_base_data", [10000], indirect=True)
def test_structure_features_regression(sampled_base_data, target_genes, gene_to_data, dataframe_regression):
    data = sampled_base_data.copy()

    with Timer("[In test] Populate DF with Structure Features"):
        processed_data = get_populated_df_with_structure_features(data, target_genes, gene_to_data)

    features = [
        STRUCTURE_SENSE_START,
        STRUCTURE_SENSE_START_FROM_END,
        STRUCTURE_SENSE_LENGTH,
        STRUCT_SENSE_IN_EXON,
        STRUCT_SENSE_IN_INTRON,
        STRUCT_SENSE_IN_UTR,
        STRUCT_SENSE_IN_CDS,
        STRUCT_SENSE_IN_CDS_NON_EXCLUSIVE,
        STRUCTURE_SENSE_TYPE,
    ]

    dataframe_regression.check(processed_data[features])
