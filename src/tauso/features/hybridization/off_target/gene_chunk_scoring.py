"""Chunk ASOs per target gene and score them across a thread pool.

The orchestration behind the per-gene off-target / on-target features: build a target
FASTA per gene, split each gene's ASOs into worker-sized chunks, score the chunks in
parallel, and assemble one score column per energy cutoff. The hit-level scoring itself
is a sum(exp(-energy/RT)) per trigger (see sum_exp_by_trigger_multi_cutoff).
"""

import logging
import os
import uuid
from collections import defaultdict

import pandas as pd

from ....data.consts import ASO_SEQUENCE
from ....util import get_antisense
from ..fast_hybridization import (
    TMP_PATH,
    Interaction,
    dump_target_file,
    parse_risearch_hits_pyarrow,
    sum_exp_by_trigger_multi_cutoff,
)
from .parallel import run_tasks_parallel

logger = logging.getLogger(__name__)

_RT = 0.616


def _validate_genes_found(target_genes, gene_to_data):
    not_found = [g for g in target_genes if g not in gene_to_data]
    if not_found:
        raise ValueError(f"The following genes are not found in gene_to_data: {not_found}")


def chunkify_score_one_gene(row_triggers, target_path, cutoffs, aggregation_factory=None, chunk_size=100):
    """Score one gene's ASOs against its target, in query chunks, for every cutoff.

    row_triggers: [(aso_index, trigger_seq)] for ASOs targeting this gene.
    ``aggregation_factory(cutoffs)`` builds the RIsearch aggregation; the default is a per-trigger
    ``sum(exp(-energy/RT))`` so the return is ``{cutoff: {trigger_id: score}}``. A different factory
    (e.g. ``stats_by_trigger_multi_cutoff``) yields tuple-valued per-trigger entries instead.
    chunk_size caps the per-call query batch (and thus peak RIsearch output held while parsing).
    """
    cutoffs = [int(c) for c in cutoffs]
    minimum_score = min(cutoffs)
    if aggregation_factory is None:
        aggregation_factory = lambda cs: sum_exp_by_trigger_multi_cutoff(cs, RT=_RT)
    combined: dict = {c: {} for c in cutoffs}
    for i in range(0, len(row_triggers), chunk_size):
        chunk = row_triggers[i : i + chunk_size]
        partial = parse_risearch_hits_pyarrow(
            trigger_id_seq_pairs=[(str(idx), trig) for idx, trig in chunk],
            target_file_path=target_path,
            aggregation=aggregation_factory(cutoffs),
            minimum_score=minimum_score,
            parsing_type="2",
            interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
            transpose=True,
            batch_id=f"{os.getpid()}-{uuid.uuid4().hex}",
        )
        # Each trigger id (an ASO row index) appears in exactly one chunk, so the per-cutoff dicts are
        # disjoint across chunks and update is exact for both scalar and tuple aggregation values.
        for c in cutoffs:
            combined[c].update(partial.get(c, {}))
    return combined


def build_gene_chunk_tasks(gene_to_row_triggers, gene_to_target_path, n_jobs):
    """Split each gene's ASOs into sub-batches so all n_jobs workers are used.

    When there are fewer genes than workers, each gene's rows are divided into
    ceil(n_jobs / n_genes) sub-batches. Returns [(gene, sub_triggers, target_path)].
    """
    n_genes = max(1, len(gene_to_row_triggers))
    tasks_per_gene = max(1, (n_jobs + n_genes - 1) // n_genes)
    tasks = []
    for gene, row_triggers in gene_to_row_triggers.items():
        sub_size = max(1, (len(row_triggers) + tasks_per_gene - 1) // tasks_per_gene)
        target_path = gene_to_target_path[gene]
        for i in range(0, len(row_triggers), sub_size):
            tasks.append((gene, row_triggers[i : i + sub_size], target_path))
    return tasks


def dispatch_gene_chunk_scores(
    aso_df,
    gene_to_data,
    target_genes,
    get_gene_fn,
    feature_name_fn,
    cutoffs,
    n_jobs,
    aggregation_factory=None,
    build_columns=None,
):
    """Score ASOs against their target gene(s) in parallel, one feature column per cutoff.

    get_gene_fn(row) returns the target gene for a row (its own canonical gene for
    on-target, or a fixed gene for single off-target). feature_name_fn(cutoff) names each
    output column. Returns (aso_df, [feature_names]).

    ``aggregation_factory`` (default: per-trigger sum-exp) selects the RIsearch aggregation.
    ``build_columns(aso_df, gene_scores, gene_to_row_triggers, cutoffs) -> [feature_names]`` overrides
    the default one-scalar-column-per-cutoff assembly, e.g. to emit several columns from a tuple-valued
    aggregation. When given, ``feature_name_fn`` is ignored.
    """
    cutoffs = [int(c) for c in cutoffs]
    _validate_genes_found(target_genes, gene_to_data)
    TMP_PATH.mkdir(exist_ok=True)

    gene_to_target_path = {}
    try:
        for gene in target_genes:
            gene_to_target_path[gene] = dump_target_file(
                f"target-{gene}-{uuid.uuid4().hex}.fa", {gene: gene_to_data[gene].full_mrna}
            )

        gene_to_row_triggers = defaultdict(list)
        for idx, seq, gene in zip(aso_df.index, aso_df[ASO_SEQUENCE], aso_df.apply(get_gene_fn, axis=1)):
            if pd.notna(gene) and gene in gene_to_target_path:
                gene_to_row_triggers[gene].append((idx, get_antisense(seq)))

        chunk_tasks = build_gene_chunk_tasks(gene_to_row_triggers, gene_to_target_path, n_jobs)
        logger.info(
            "RIsearch scoring: %d genes, %d ASOs, %d tasks, %d cutoffs",
            len(gene_to_row_triggers),
            sum(len(v) for v in gene_to_row_triggers.values()),
            len(chunk_tasks),
            len(cutoffs),
        )

        tasks = [
            ((gene, i), triggers, path, cutoffs, aggregation_factory)
            for i, (gene, triggers, path) in enumerate(chunk_tasks)
        ]
        chunk_scores = run_tasks_parallel(tasks, chunkify_score_one_gene, n_jobs)

        gene_scores: dict = defaultdict(lambda: {c: {} for c in cutoffs})
        for (gene, _), per_cutoff in chunk_scores.items():
            for c in cutoffs:
                gene_scores[gene][c].update(per_cutoff.get(c, {}))

        if build_columns is not None:
            feature_names = build_columns(aso_df, gene_scores, gene_to_row_triggers, cutoffs)
        else:
            feature_names = []
            for cutoff in cutoffs:
                scores = pd.Series(0.0, index=aso_df.index)
                for gene, row_triggers in gene_to_row_triggers.items():
                    gene_cutoff_scores = gene_scores[gene][cutoff]
                    for idx, _ in row_triggers:
                        value = gene_cutoff_scores.get(str(idx))
                        if value is not None:
                            scores[idx] = value
                fname = feature_name_fn(cutoff)
                aso_df[fname] = scores
                feature_names.append(fname)
    finally:
        for path in gene_to_target_path.values():
            if os.path.exists(path):
                os.remove(path)

    return aso_df, feature_names
