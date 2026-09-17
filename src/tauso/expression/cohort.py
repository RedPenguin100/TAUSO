"""Per-cell-line expression, as the features want it: {depmap_id: DataFrame}."""

import os

from ..common.gtf import filter_gtf_genes
from ..data.data import get_data_dir, load_gtf_db
from ..dependencies.depmap import (
    load_cell_line_gene_expression,
    load_cell_line_gene_transcripts,
    load_cell_line_transcript_expression,
)

GENE_EXPRESSION_DIR = "processed_expression"
TRANSCRIPT_EXPRESSION_DIR = "processed_transcript_expression"


def _dir(name):
    return os.path.join(get_data_dir(), name)


def load_cohort_expression(cell_lines_depmap, genome="GRCh38"):
    """Gene-level expression per cell line, for the genes the annotation calls valid.

    Mitochondrial genes are left out: they are expressed far above everything else and are
    not what an ASO is designed against.
    """
    valid_genes = filter_gtf_genes(load_gtf_db(genome), filter_mode="non_mt")
    return load_cell_line_gene_expression(cell_lines_depmap, valid_genes, expression_dir=_dir(GENE_EXPRESSION_DIR))


def load_cohort_transcript_expression(cell_lines_depmap, transcript_names):
    """Transcript-level expression per cell line, for the named transcripts only.

    The DepMap table has ~237,000 transcripts and the features name a handful of them.
    """
    return load_cell_line_transcript_expression(
        cell_lines_depmap, transcript_names, expression_dir=_dir(TRANSCRIPT_EXPRESSION_DIR)
    )


def load_cohort_gene_transcripts(cell_lines_depmap, genes):
    """Transcript-level expression per cell line, for whole genes rather than named
    transcripts: how a gene's expression is spread across its isoforms."""
    return load_cell_line_gene_transcripts(cell_lines_depmap, genes, expression_dir=_dir(TRANSCRIPT_EXPRESSION_DIR))
