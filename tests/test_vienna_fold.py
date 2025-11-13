import pytest
import pickle
import numpy as np

from tauso.fold import get_weighted_energy, calculate_energies
from tests.conftest import TEST_CACHE_PATH
from tauso.read_human_genome import get_locus_to_data_dict


@pytest.fixture
def mrna():
    target_gene = 'DDX11L1'

    test_cache = TEST_CACHE_PATH / 'gene_to_data_test.pickle'
    if not test_cache.exists():
        gene_to_data = get_locus_to_data_dict(include_introns=True, gene_subset=[target_gene])
        with open(test_cache, 'wb') as f:
            pickle.dump(gene_to_data, f)
    else:
        with open(test_cache, 'rb') as f:
            gene_to_data = pickle.load(f)

    yield gene_to_data[target_gene].full_mrna


def test_regression(mrna):
    energies = calculate_energies(str(mrna), 15, 40)
    print(energies[-3:])
    energy = get_weighted_energy(2525, 16, 15, energies, 40)

    assert pytest.approx(energy, rel=1e-2) == -3.56



def test_sanity(mrna):
    gene_length = 2541

    window_size = 40
    step_size = 15
    sample_energies = np.zeros((gene_length - window_size) // step_size + 1 + 1)
    sample_energies[-1:] = -1


    energy = get_weighted_energy(2525, 16, 15, sample_energies, 40)
    # Last window ends in 2540, unique coverage for 2530-2541
    # Before last window is 0, but joint with last.
    # So in total the fold should be the weighted average of the two windows
    assert pytest.approx(energy, rel=1e-2) == (5 * ((0 - 1) / 2) + 11 * (-1) ) / 16

    energy = get_weighted_energy(2526, 16, 15, sample_energies, 40)
    assert pytest.approx(energy, rel=1e-2) == (4 * ((0 - 1) / 2) + 12 * (-1) ) / 16

    sample_energies[-2:-1] = -2
    energy = get_weighted_energy(2525, 16, 15, sample_energies, 40)
    assert pytest.approx(energy, rel=1e-2) == (5 * ((-2 - 1) / 2) + 11 * (-1) ) / 16

    sample_energies[-2:] = -2
    energy = get_weighted_energy(2525, 16, 15, sample_energies, 40)
    assert pytest.approx(energy, rel=1e-2) == -2