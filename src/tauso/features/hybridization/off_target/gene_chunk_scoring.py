"""Per-gene RIsearch scan engine.

Scores ASOs against their target gene(s) with one loose RIsearch pass per gene, in parallel: a
target FASTA per gene, ASOs split into worker chunks over a thread pool, hits reduced to
a table of ``EnergyStats`` fields per ASO and cutoff. That superset feeds every
per-gene hybridization feature, each of which reads the columns it needs
(see off_target_specific_gene.py).
"""

import logging
from collections import defaultdict
from math import ceil
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd
import pyrisearch_tauso
from pyrisearch_tauso import EnergyStats

from ....data.consts import ASO_SEQUENCE
from . import ASO_TARGET_MATRIX, EXTENSION_PENALTY, RT_KCAL_MOL
from .parallel import run_tasks_parallel

logger = logging.getLogger(__name__)


def _validate_genes_found(target_genes, gene_to_data):
    not_found = [g for g in target_genes if g not in gene_to_data]
    if not_found:
        raise ValueError(f"The following genes are not found in gene_to_data: {not_found}")


def _scan_one_gene_task(row_queries, target_path, cutoffs):
    """Site-resolved stats for one gene's ASOs against its target, per cutoff.

    row_queries: [(aso_index, query_seq)] for ASOs targeting this gene.
    Returns ``{cutoff: {query_id: EnergyStats}}``.
    """
    if not row_queries:
        return {c: {} for c in cutoffs}
    return pyrisearch_tauso.search_reduced(
        queries=[(str(idx), seq) for idx, seq in row_queries],
        targets=target_path,
        reduction=pyrisearch_tauso.energy_stats(cutoffs, group_by=("query",), rt=RT_KCAL_MOL),
        min_score=min(cutoffs),
        matrix=ASO_TARGET_MATRIX,
        extension_penalty=EXTENSION_PENALTY,
        transpose=True,
    )


def build_gene_chunk_tasks(gene_to_row_queries, gene_to_target_path, n_jobs):
    """One search per gene, split further so that every worker has something to do.

    A gene is one search, so a run over few genes would leave most workers idle;
    each gene's ASOs are divided into enough batches to go round. Returns
    [(queries, target_path)].
    """
    batches_per_gene = ceil(n_jobs / max(1, len(gene_to_row_queries)))
    tasks = []
    for gene, row_queries in gene_to_row_queries.items():
        batch = max(1, ceil(len(row_queries) / batches_per_gene))
        for start in range(0, len(row_queries), batch):
            tasks.append((row_queries[start : start + batch], gene_to_target_path[gene]))
    return tasks


def scan_gene_sites(aso_df, gene_to_data, target_genes, row_genes, cutoffs, n_jobs):
    """Search each ASO against the gene it is assigned, and tabulate what was found.

    row_genes gives the target gene of each row, in DataFrame order. The result has
    one row per ASO and a column per (cutoff, statistic) -- sum_exp, min_energy and
    n_sites, the fields of EnergyStats. An ASO with no hit above a cutoff is NaN
    there, so each feature decides for itself what that means.
    """
    cutoffs = [int(c) for c in cutoffs]
    _validate_genes_found(target_genes, gene_to_data)

    with TemporaryDirectory() as work:
        target_path = {}
        for i, gene in enumerate(target_genes):
            path = Path(work) / f"{i}.fa"
            path.write_text(f">{gene}\n{gene_to_data[gene].full_mrna}\n")
            target_path[gene] = str(path)

        gene_to_row_queries = defaultdict(list)
        for idx, seq, gene in zip(aso_df.index, aso_df[ASO_SEQUENCE], row_genes):
            if pd.notna(gene) and gene in target_path:
                gene_to_row_queries[gene].append((idx, seq))

        tasks = build_gene_chunk_tasks(gene_to_row_queries, target_path, n_jobs)
        logger.info(
            "RIsearch scan: %d genes, %d ASOs, %d tasks, %d cutoffs",
            len(gene_to_row_queries),
            sum(len(v) for v in gene_to_row_queries.values()),
            len(tasks),
            len(cutoffs),
        )
        results = run_tasks_parallel(
            [(i, queries, path, cutoffs) for i, (queries, path) in enumerate(tasks)],
            _scan_one_gene_task,
            n_jobs,
        )

    # A query is named for the row it came from; this reads the name back.
    row_of = {str(idx): idx for idx in aso_df.index}
    columns = {(cutoff, field): {} for cutoff in cutoffs for field in EnergyStats._fields}
    for per_cutoff in results.values():
        for cutoff, by_query in per_cutoff.items():
            for query_id, stats in by_query.items():
                row = row_of[query_id]
                for field, value in zip(EnergyStats._fields, stats):
                    columns[cutoff, field][row] = value

    return pd.DataFrame(columns, index=aso_df.index, dtype=float)
