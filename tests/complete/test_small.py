import pytest

from tauso.features.names import *
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
        SENSE_START,
        SENSE_START_FROM_END,
        SENSE_LENGTH,
        SENSE_EXON,
        SENSE_INTRON,
        SENSE_UTR,
        SENSE_TYPE,
        "n_sense_3utr",
        "n_sense_5utr",
        "n_sense_exon",
        "n_sense_intron",
        "n_sense_utr",
    ]

    dataframe_regression.check(processed_data[features])
