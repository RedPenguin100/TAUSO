import logging
import os
import uuid
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

logger = logging.getLogger(__name__)

from ....data.consts import CANONICAL_GENE, SEQUENCE
from ....util import get_antisense
from ..fast_hybridization import (
    TMP_PATH,
    Interaction,
    dump_target_file,
    parse_risearch_hits_pyarrow,
    sum_exp_by_trigger_multi_cutoff,
)

_RT = 0.616


def _validate_genes_found(target_genes, gene_to_data):
    not_found = [g for g in target_genes if g not in gene_to_data]
    if not_found:
        raise ValueError(f"The following genes are not found in gene_to_data: {not_found}")


def _score_one_gene_multi_cutoff(gene, row_triggers, target_path, cutoffs, chunk_size=100):
    """Run RIsearch for one gene in chunks at the loosest cutoff and derive the
    sum(exp(-RT·energy)) score for every cutoff from that single pass. Thread-safe.

    Returns {cutoff: {trigger_id: score}}. Uses the pyarrow streaming parser
    (sum_exp_by_trigger_multi_cutoff): the parse + exp + group_by all run in C++ and
    release the GIL, so this scales under the ThreadPoolExecutor in _apply_risearch_scoring.
    Memory is bounded by the pyarrow read block, not the (multi-hundred-MB at cutoff=0)
    output. chunk_size caps peak memory and the per-call query batch.
    """
    cutoffs = [int(c) for c in cutoffs]
    combined: dict = {c: {} for c in cutoffs}
    minimum_score = min(cutoffs)
    for i in range(0, len(row_triggers), chunk_size):
        chunk = row_triggers[i : i + chunk_size]
        partial = parse_risearch_hits_pyarrow(
            trigger_id_seq_pairs=[(str(idx), trig) for idx, trig in chunk],
            target_file_path=target_path,
            aggregation=sum_exp_by_trigger_multi_cutoff(cutoffs, rt=_RT),
            minimum_score=minimum_score,
            parsing_type="2",
            interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
            transpose=True,
            batch_id=f"{os.getpid()}-{uuid.uuid4().hex}",
        )
        # Sum across chunks for the same trigger (rare, only if a trigger straddles
        # the chunksize boundary), per cutoff.
        for c in cutoffs:
            dst = combined[c]
            for trig, v in partial.get(c, {}).items():
                dst[trig] = dst.get(trig, 0.0) + v
    return combined


def _build_tasks(gene_to_row_triggers, gene_to_target_info, n_jobs):
    """
    Split each gene's row-triggers into sub-lists so that all n_jobs workers are used.

    When n_jobs exceeds the number of genes, each gene's rows are divided into
    ceil(n_jobs / n_genes) sub-tasks so the spare workers aren't wasted.
    Returns a list of (gene, sub_triggers, target_path) tuples.
    """
    n_genes = max(1, len(gene_to_row_triggers))
    tasks_per_gene = max(1, (n_jobs + n_genes - 1) // n_genes)
    tasks = []
    for gene, row_triggers in gene_to_row_triggers.items():
        sub_size = max(1, (len(row_triggers) + tasks_per_gene - 1) // tasks_per_gene)
        target_path = gene_to_target_info[gene]
        for i in range(0, len(row_triggers), sub_size):
            tasks.append((gene, row_triggers[i : i + sub_size], target_path))
    return tasks


def _apply_risearch_scoring(
    aso_df,
    gene_to_data,
    target_genes,
    get_gene_fn,
    feature_name_fn,
    cutoffs,
    n_jobs,
    verbose,
    stream=True,
):
    """Core logic for RIsearch hybridization scoring, producing one feature per cutoff.

    Each (gene, ASO sub-batch) tile runs RIsearch ONCE at the loosest cutoff and derives
    every cutoff's score from that single streaming pass (cutoff-collapse). Tiles run on
    a ThreadPoolExecutor. `feature_name_fn(cutoff)` names each output column; `stream` is
    accepted for backward compatibility (the pyarrow multi-cutoff path always streams).
    """
    cutoffs = [int(c) for c in cutoffs]
    _validate_genes_found(target_genes, gene_to_data)
    TMP_PATH.mkdir(exist_ok=True)

    gene_to_target_info = {}
    try:
        for gene in target_genes:
            sequence = gene_to_data[gene].full_mrna
            target_path = dump_target_file(f"target-{gene}-{uuid.uuid4().hex}.fa", {gene: sequence})
            gene_to_target_info[gene] = target_path

        gene_to_row_triggers = defaultdict(list)
        for idx, seq, gene in zip(
            aso_df.index,
            aso_df[SEQUENCE],
            aso_df.apply(get_gene_fn, axis=1),
        ):
            if pd.isna(gene) or gene not in gene_to_target_info:
                continue
            gene_to_row_triggers[gene].append((idx, get_antisense(seq)))

        tasks = _build_tasks(gene_to_row_triggers, gene_to_target_info, n_jobs)
        effective_workers = min(n_jobs, len(tasks)) if n_jobs > 1 else 1
        logger.info(
            "RIsearch scoring: %d genes × %d rows → %d tasks across %d workers, %d cutoffs",
            len(gene_to_row_triggers),
            sum(len(v) for v in gene_to_row_triggers.values()),
            len(tasks),
            effective_workers,
            len(cutoffs),
        )

        # Accumulate per gene, per cutoff: gene_scores[gene][cutoff] = {trigger: score}
        gene_scores: dict = defaultdict(lambda: {c: {} for c in cutoffs})

        def _merge(gene, per_cutoff):
            for c in cutoffs:
                gene_scores[gene][c].update(per_cutoff.get(c, {}))

        if effective_workers > 1:
            with ThreadPoolExecutor(max_workers=effective_workers) as pool:
                futures = [
                    (pool.submit(_score_one_gene_multi_cutoff, gene, sub_triggers, target_path, cutoffs), gene)
                    for gene, sub_triggers, target_path in tasks
                ]
                for fut, gene in futures:
                    _merge(gene, fut.result())
        else:
            for gene, sub_triggers, target_path in tasks:
                _merge(gene, _score_one_gene_multi_cutoff(gene, sub_triggers, target_path, cutoffs))

        feature_names = []
        for cutoff in cutoffs:
            scores = pd.Series(0.0, index=aso_df.index)
            for gene, row_triggers in gene_to_row_triggers.items():
                score_dict = gene_scores[gene][cutoff]
                for idx, _ in row_triggers:
                    val = score_dict.get(str(idx))
                    if val is not None:
                        scores[idx] = val
            fname = feature_name_fn(cutoff)
            aso_df[fname] = scores
            feature_names.append(fname)

    finally:
        for path in gene_to_target_info.values():
            if os.path.exists(path):
                os.remove(path)

    return aso_df, feature_names


def on_target_total_hybridization(aso_df, gene_to_data, cutoffs, n_jobs=1, verbose=False, stream=True):
    """Scores each oligo against its own canonical gene, one feature per cutoff.

    All cutoffs are derived from a single loose RIsearch pass per (gene, ASO batch); the
    pyarrow streaming parser keeps memory bounded regardless of cutoff (notably cutoff=0,
    whose materialized TSV can be hundreds of MB). Returns (df, [feature_names]).
    """
    unique_genes = aso_df[CANONICAL_GENE].dropna().unique()
    return _apply_risearch_scoring(
        aso_df=aso_df,
        gene_to_data=gene_to_data,
        target_genes=unique_genes,
        get_gene_fn=lambda row: row[CANONICAL_GENE],
        feature_name_fn=lambda c: f"on_target_total_hybridization_{c}",
        cutoffs=cutoffs,
        n_jobs=n_jobs,
        verbose=verbose,
        stream=stream,
    )


def off_target_specific_seq_pandarallel(aso_df, gene_name, gene_to_data, cutoffs, n_jobs=1, verbose=False, stream=True):
    """Scores each oligo against a single statically provided gene, one feature per cutoff.

    All cutoffs are derived from a single loose RIsearch pass per ASO batch (cutoff-collapse),
    on a ThreadPoolExecutor + pyarrow streaming parser. Returns (df, [feature_names]).
    (The ``_pandarallel`` name is legacy — this path has not used pandarallel for a while.)
    """
    return _apply_risearch_scoring(
        aso_df=aso_df,
        gene_to_data=gene_to_data,
        target_genes=[gene_name],
        get_gene_fn=lambda row: gene_name,
        feature_name_fn=lambda c: f"off_target_single_{gene_name}_c{c}",
        cutoffs=cutoffs,
        n_jobs=n_jobs,
        verbose=verbose,
        stream=stream,
    )
