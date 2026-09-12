"""Per-gene RIsearch scan engine.

Scores ASOs against their target gene(s) with one loose RIsearch pass per gene, in parallel: a
target FASTA per gene, ASOs split into worker chunks over a thread pool, hits reduced to
site-resolved ``EnergyStats`` (Boltzmann sum, best-site energy) per ASO and cutoff.
That superset feeds every per-gene hybridization feature, which derive their columns via
``emit_site_columns`` (see off_target_specific_gene.py).
"""

import logging
from collections import defaultdict
from contextlib import ExitStack

import pandas as pd
import pyrisearch_tauso

from ....data.consts import ASO_SEQUENCE
from ....pandas_utils import add_columns
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
    cutoffs = [int(c) for c in cutoffs]
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
    """Split each gene's ASOs into sub-batches so all n_jobs workers are used.

    When there are fewer genes than workers, each gene's rows are divided into
    ceil(n_jobs / n_genes) sub-batches. Returns [(gene, sub_queries, target_path)].
    """
    n_genes = max(1, len(gene_to_row_queries))
    tasks_per_gene = max(1, (n_jobs + n_genes - 1) // n_genes)
    tasks = []
    for gene, row_queries in gene_to_row_queries.items():
        sub_size = max(1, (len(row_queries) + tasks_per_gene - 1) // tasks_per_gene)
        target_path = gene_to_target_path[gene]
        for i in range(0, len(row_queries), sub_size):
            tasks.append((gene, row_queries[i : i + sub_size], target_path))
    return tasks


def scan_gene_sites(aso_df, gene_to_data, target_genes, row_genes, cutoffs, n_jobs):
    """Run the per-gene RIsearch scan and return site-resolved stats per ASO and cutoff.

    row_genes supplies the target gene for each row, in DataFrame order.
    Returns ``{cutoff: {row_index: EnergyStats}}`` covering every ASO that had
    at least one hit above the cutoff against its target.
    """
    cutoffs = [int(c) for c in cutoffs]
    _validate_genes_found(target_genes, gene_to_data)

    with ExitStack() as targets:
        gene_to_target_path = {
            gene: targets.enter_context(pyrisearch_tauso.fasta_targets({gene: gene_to_data[gene].full_mrna}))
            for gene in target_genes
        }

        gene_to_row_queries = defaultdict(list)
        for idx, seq, gene in zip(aso_df.index, aso_df[ASO_SEQUENCE], row_genes):
            if pd.notna(gene) and gene in gene_to_target_path:
                gene_to_row_queries[gene].append((idx, seq))

        chunk_tasks = build_gene_chunk_tasks(gene_to_row_queries, gene_to_target_path, n_jobs)
        logger.info(
            "RIsearch scan: %d genes, %d ASOs, %d tasks, %d cutoffs",
            len(gene_to_row_queries),
            sum(len(v) for v in gene_to_row_queries.values()),
            len(chunk_tasks),
            len(cutoffs),
        )

        tasks = [((gene, i), queries, path, cutoffs) for i, (gene, queries, path) in enumerate(chunk_tasks)]
        chunk_stats = run_tasks_parallel(tasks, _scan_one_gene_task, n_jobs)

        scan = {c: {} for c in cutoffs}
        for (_gene, i), per_cutoff in chunk_stats.items():
            for idx, _ in chunk_tasks[i][1]:
                for c in cutoffs:
                    stat = per_cutoff.get(c, {}).get(str(idx))
                    if stat is not None:
                        scan[c][idx] = stat

    return scan


def emit_site_columns(aso_df, scan, cutoffs, derivations):
    """Write one column per (derivation, cutoff) from a scan_gene_sites result.

    derivations: [(name_fn, derive_fn)]. For each cutoff, name_fn(cutoff) names the column and
    derive_fn(EnergyStats) -> value; rows with no scored hit (or a None derivation) default to 0.0.
    Returns (aso_df, [feature_names]).
    """
    cutoffs = [int(c) for c in cutoffs]
    columns = {}
    for cutoff in cutoffs:
        per_row = scan.get(cutoff, {})
        for name_fn, derive_fn in derivations:
            derived = {idx: value for idx, stats in per_row.items() if (value := derive_fn(stats)) is not None}
            mapped = pd.Series(aso_df.index.map(derived.get), index=aso_df.index, dtype=float)
            columns[name_fn(cutoff)] = mapped.fillna(0.0)

    return add_columns(aso_df, columns), list(columns)
