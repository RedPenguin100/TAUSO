import os

from ..common.gtf import filter_gtf_genes
from ..data.data import get_data_dir, load_gtf_db
from ..dependencies.depmap import load_cell_line_gene_expression


def load_transcriptomes(cell_lines_depmap):
    # 1. Load the full database and get valid non-mitochondrial genes
    db = load_gtf_db()
    valid_genes = filter_gtf_genes(db, filter_mode="non_mt")

    data_dir = get_data_dir()
    expression_dir = os.path.join(data_dir, "processed_expression")

    # 2. Load cell-line specific transcriptomes

    transcriptomes = load_cell_line_gene_expression(
        cell_lines_depmap,
        valid_genes,
        expression_dir=expression_dir,
    )

    return transcriptomes
