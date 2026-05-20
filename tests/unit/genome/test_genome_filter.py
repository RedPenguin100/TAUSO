import logging

import pytest

from tauso.common.gtf import filter_gff_genes, filter_gtf_genes
from tauso.data.data import load_gff_db, load_gtf_db

logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


@pytest.fixture(params=["gtf", "gff"])
def db_and_filter(request):
    if request.param == "gtf":
        return load_gtf_db(), filter_gtf_genes
    else:
        return load_gff_db(), filter_gff_genes


@pytest.mark.parametrize("filter_mode", ["non_mt", "protein_coding"])
def test_filter_genes(db_and_filter, filter_mode, data_regression):
    db, filter_fn = db_and_filter
    result = filter_fn(db, filter_mode=filter_mode)
    data_regression.check(result)


@pytest.mark.parametrize("filter_mode", ["non_mt", "protein_coding"])
def test_filter_genes_gtf_gff_match(filter_mode):
    gtf_result = filter_gtf_genes(load_gtf_db(), filter_mode=filter_mode)
    gff_result = filter_gff_genes(load_gff_db(), filter_mode=filter_mode)
    assert gtf_result == gff_result, (
        f"GTF/GFF mismatch for {filter_mode}: GTF-only={gtf_result - gff_result}, GFF-only={gff_result - gtf_result}"
    )
