import os
from pathlib import Path

from ..common.gtf import filter_gtf_genes
from ..data.data import get_data_dir, load_db
from ..features.codon_usage.find_cai_reference import load_cell_line_gene_expression
from ..features.hybridization_off_target.common import get_general_expression_of_genes


def load_transcriptomes(cell_lines_depmap):
    # 1. Load the full database and get valid non-mitochondrial genes
    db = load_db()
    valid_genes = filter_gtf_genes(db, filter_mode="non_mt")

    data_dir = get_data_dir()
    expression_dir = os.path.join(data_dir, "processed_expression")

    # 2. Load cell-line specific transcriptomes

    transcriptomes = load_cell_line_gene_expression(
        cell_lines_depmap,
        valid_genes,
        expression_dir=expression_dir,
    )

    # 3. Load the general mean expression
    mean_exp_data = get_general_expression_of_genes(
        Path(data_dir) / "OmicsExpressionTPMLogp1HumanAllGenesStranded.csv", valid_genes
    )
    transcriptomes["general"] = mean_exp_data
    return transcriptomes
