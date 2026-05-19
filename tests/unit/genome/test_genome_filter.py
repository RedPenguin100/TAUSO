import logging

from tauso.common.gtf import filter_gtf_genes
from tauso.data.data import load_db

# Configure logging to show debug messages in your notebook
logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def test_filter_genes_non_mt(data_regression):
    pr = load_db()
    logger.debug("Loaded GTF")
    non_mt = filter_gtf_genes(pr, filter_mode="non_mt")
    data_regression.check(non_mt)


def test_filter_genes_protein_coding(data_regression):
    pr = load_db()
    logger.debug("Loaded GTF")
    protein_coding = filter_gtf_genes(pr, filter_mode="protein_coding")
    data_regression.check(protein_coding)
