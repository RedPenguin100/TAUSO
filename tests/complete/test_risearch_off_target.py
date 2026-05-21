import pytest
from notebooks.data.OligoAI.parse_chemistry import assign_chemistry

from tests.complete.conftest import get_n_jobs
from tauso.features.hybridization_off_target.off_target_feature import populate_off_target_specific
from tauso.genome.transcriptome import load_transcriptomes
from tauso.data.consts import CELL_LINE_DEPMAP


@pytest.fixture(scope="module")
def chemistry_data(structure_data):
    return assign_chemistry(structure_data)


@pytest.fixture
def mini_sampled_data(request, chemistry_data):
    """Samples the fully processed DataFrame right before the test runs."""
    n_samples = getattr(request, "param", 1000)
    actual_samples = min(n_samples, len(chemistry_data))
    return chemistry_data.sample(n=actual_samples, random_state=42).copy()


@pytest.mark.parametrize("mini_sampled_data", [10], indirect=True)
def test_off_target_hybridization(request, mini_sampled_data, dataframe_regression):
    gene_to_data_full = request.getfixturevalue("gene_to_data_full")

    cell_lines_depmap = mini_sampled_data[CELL_LINE_DEPMAP].dropna().unique().tolist()
    transcriptomes = load_transcriptomes(cell_lines_depmap)

    mini_sampled_data, feature_names = populate_off_target_specific(
        mini_sampled_data, gene_to_data_full, transcriptomes, [50], [1200], method="ARTM", n_jobs=get_n_jobs()
    )

    dataframe_regression.check(mini_sampled_data[["index_oligo"] + feature_names])
