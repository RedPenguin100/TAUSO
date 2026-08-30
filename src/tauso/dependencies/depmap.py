import logging
import os

import pandas as pd

logger = logging.getLogger(__name__)


def load_cell_line_gene_expression(depmap_ids, valid_genes, expression_dir):
    transcriptomes = {}

    for ach_id in depmap_ids:
        path = os.path.join(expression_dir, f"{ach_id}_expression.csv")
        if not os.path.exists(path):
            logger.warning("No expression data for %s", ach_id)
            continue

        # Load and Filter
        df = pd.read_csv(path)
        df = df[df["Gene"].isin(valid_genes)].copy()

        if not df.empty:
            transcriptomes[ach_id] = df

    return transcriptomes


def load_cell_line_transcript_expression(depmap_ids, valid_transcript_names, expression_dir):
    """Per-cell-line transcript expression, the transcript twin of
    load_cell_line_gene_expression.

    Reads the files build-cohort-transcript-expression writes, filtered to the Ensembl
    transcript names asked for (RNASEH1-201 and the like), which say which gene and which
    isoform without a lookup.
    """
    wanted = set(valid_transcript_names)
    transcriptomes = {}

    for ach_id in depmap_ids:
        path = os.path.join(expression_dir, f"{ach_id}_transcript_expression.csv")
        if not os.path.exists(path):
            logger.warning("No transcript expression data for %s", ach_id)
            continue

        df = pd.read_csv(path)
        df = df[df["TranscriptName"].isin(wanted)].copy()

        if not df.empty:
            transcriptomes[ach_id] = df

    return transcriptomes


def load_cell_line_gene_transcripts(depmap_ids, valid_genes, expression_dir):
    """Per-cell-line transcript expression for whole genes rather than named transcripts.

    Where load_cell_line_transcript_expression answers "how much of this one isoform", this
    answers "how is this gene's expression distributed across its isoforms", which is what a
    feature about the target gene needs: the target differs from row to row, so the isoforms
    cannot be named ahead of time.
    """
    wanted = set(valid_genes)
    transcriptomes = {}

    for ach_id in depmap_ids:
        path = os.path.join(expression_dir, f"{ach_id}_transcript_expression.csv")
        if not os.path.exists(path):
            logger.warning("No transcript expression data for %s", ach_id)
            continue

        df = pd.read_csv(path)
        df = df[df["Gene"].isin(wanted)].copy()

        if not df.empty:
            transcriptomes[ach_id] = df

    return transcriptomes
