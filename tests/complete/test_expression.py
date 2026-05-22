import os

import pytest

from tauso.data.consts import CELL_LINE_DEPMAP
from tauso.data.data import get_data_dir
from tauso.features.codon_usage.find_cai_reference import load_cell_line_gene_expression
from tauso.populate.populate_context import (
    _SPECIAL_GENES,
    populate_special_gene_expression,
    populate_target_expression,
)
from tauso.timer import Timer


@pytest.fixture(scope="session")
def expression_transcriptomes(base_data, target_genes):
    """Expression data for target genes plus the fixed special genes (RNASEH1, STAB2, MRC1, MSR1)."""
    data_dir = get_data_dir()
    expression_dir = os.path.join(data_dir, "processed_expression")
    depmap_ids = list(set(base_data[CELL_LINE_DEPMAP]))
    valid_genes = list(set(target_genes) | set(_SPECIAL_GENES.values()))
    with Timer("Load Expression Transcriptomes"):
        return load_cell_line_gene_expression(depmap_ids, valid_genes, expression_dir=expression_dir)


@pytest.mark.parametrize("mini_sampled_data", [1000], indirect=True)
def test_target_expression_regression(mini_sampled_data, expression_transcriptomes, dataframe_regression):
    data = mini_sampled_data.copy()
    result, feats = populate_target_expression(data, expression_transcriptomes)
    dataframe_regression.check(result[feats])


@pytest.mark.parametrize("mini_sampled_data", [1000], indirect=True)
def test_special_gene_expression_regression(mini_sampled_data, expression_transcriptomes, dataframe_regression):
    data = mini_sampled_data.copy()
    result, feats = populate_special_gene_expression(data, expression_transcriptomes)
    dataframe_regression.check(result[feats])
