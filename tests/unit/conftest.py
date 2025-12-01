import pytest

import pickle
from tauso.genome.read_human_genome import get_locus_to_data_dict
from tests.common.consts import UNIT_TESTS_PATH


def pytest_addoption(parser):
    parser.addoption(
        "--runslow", action="store_true", default=False, help="run slow tests"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        # --runslow given in cli: do not skip slow tests
        return
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)


@pytest.fixture
def short_mrna():
    target_gene = 'DDX11L1'

    test_cache = UNIT_TESTS_PATH / 'gene_to_data_test.pickle'
    if not test_cache.exists():
        UNIT_TESTS_PATH.mkdir(parents=True, exist_ok=True)
        gene_to_data = get_locus_to_data_dict(include_introns=True, gene_subset=[target_gene])
        with open(test_cache, 'wb') as f:
            pickle.dump(gene_to_data, f)
    else:
        with open(test_cache, 'rb') as f:
            gene_to_data = pickle.load(f)

    yield gene_to_data[target_gene]

