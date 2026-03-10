import pickle

import pytest

from tauso.genome.read_human_genome import get_locus_to_data_dict

from .common.consts import TESTS_CACHE

SHORT_GENE = "DDX11L1"
TESTED_GENES = [SHORT_GENE]


@pytest.fixture
def small_gene_to_data():
    test_cache = TESTS_CACHE / "gene_to_data_test.pickle"
    if not test_cache.exists():
        TESTS_CACHE.mkdir(parents=True, exist_ok=True)
        gene_to_data = get_locus_to_data_dict(
            include_introns=True, gene_subset=TESTED_GENES
        )
        with open(test_cache, "wb") as f:
            pickle.dump(gene_to_data, f)
    else:
        with open(test_cache, "rb") as f:
            gene_to_data = pickle.load(f)
    yield gene_to_data


@pytest.fixture
def short_mrna(small_gene_to_data):
    yield small_gene_to_data[SHORT_GENE]


def pytest_addoption(parser):
    parser.addoption(
        "--run-integration",
        action="store_true",
        default=False,
        help="Run integration tests",
    )


def pytest_runtest_setup(item):
    if "integration" in item.keywords and not item.config.getoption(
        "--run-integration"
    ):
        pytest.skip("Skipping integration tests (use --run-integration to run)")
