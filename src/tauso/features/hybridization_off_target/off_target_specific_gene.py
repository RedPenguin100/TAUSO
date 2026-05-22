import logging
import os
import uuid
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

from ...data.consts import CANONICAL_GENE, SEQUENCE
from ...util import get_antisense
from ..hybridization.fast_hybridization import (
    TMP_PATH,
    Interaction,
    dump_target_file,
    get_triggers_mfe_scores_batch,
)
from .off_target_functions import parse_risearch_output

_RT = 0.616


def _sum_exp_energy_by_trigger(result_df: pd.DataFrame) -> dict:
    """Return {trigger_id: sum(exp(-RT * energy))} for every trigger in result_df."""
    exp_vals = np.exp(-_RT * result_df["energy"].to_numpy())
    return result_df.assign(_exp=exp_vals).groupby("trigger", sort=False)["_exp"].sum().to_dict()


def _validate_genes_found(target_genes, gene_to_data):
    not_found = [g for g in target_genes if g not in gene_to_data]
    if not_found:
        raise ValueError(f"The following genes are not found in gene_to_data: {not_found}")


def _score_one_gene(gene, row_triggers, target_path, cutoff, chunk_size=100):
    """Run RIsearch for one gene in chunks and return {trigger_id: score}. Thread-safe.

    Chunking caps peak memory: with cutoff=0 a single 1000-query call can
    generate hundreds of MB of RIsearch output held in RAM before parsing.
    chunk_size=100 gives 10 calls instead of 1 while still being ~100x
    faster than the old one-subprocess-per-row approach.
    """
    combined: dict = {}
    for i in range(0, len(row_triggers), chunk_size):
        chunk = row_triggers[i : i + chunk_size]
        result = get_triggers_mfe_scores_batch(
            trigger_id_seq_pairs=[(str(idx), trig) for idx, trig in chunk],
            target_file_path=target_path,
            minimum_score=cutoff,
            parsing_type="2",
            interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
            transpose=True,
            batch_id=f"{os.getpid()}-{uuid.uuid4().hex}",
        )
        if not result.strip():
            continue
        result_df = parse_risearch_output(result)
        del result
        if result_df.empty or "energy" not in result_df.columns:
            continue
        combined.update(_sum_exp_energy_by_trigger(result_df))
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
    feature_name,
    cutoff,
    n_jobs,
    verbose,
):
    """Core logic for RIsearch hybridization scoring."""
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
            "RIsearch scoring: %d genes × %d rows → %d tasks across %d workers",
            len(gene_to_row_triggers),
            sum(len(v) for v in gene_to_row_triggers.values()),
            len(tasks),
            effective_workers,
        )

        # Accumulate partial results per gene before building the score series
        gene_scores: dict[str, dict] = defaultdict(dict)

        if effective_workers > 1:
            with ThreadPoolExecutor(max_workers=effective_workers) as pool:
                futures = [
                    (pool.submit(_score_one_gene, gene, sub_triggers, target_path, cutoff), gene)
                    for gene, sub_triggers, target_path in tasks
                ]
                for fut, gene in futures:
                    gene_scores[gene].update(fut.result())
        else:
            for gene, sub_triggers, target_path in tasks:
                gene_scores[gene].update(_score_one_gene(gene, sub_triggers, target_path, cutoff))

        scores = pd.Series(0.0, index=aso_df.index)
        for gene, row_triggers in gene_to_row_triggers.items():
            score_dict = gene_scores[gene]
            for idx, _ in row_triggers:
                val = score_dict.get(str(idx))
                if val is not None:
                    scores[idx] = val

    finally:
        for path in gene_to_target_info.values():
            if os.path.exists(path):
                os.remove(path)

    aso_df[feature_name] = scores
    return aso_df, feature_name


def on_target_total_hybridization(aso_df, gene_to_data, cutoff, n_jobs=1, verbose=False):
    """Scores each oligo against its own canonical gene."""
    unique_genes = aso_df[CANONICAL_GENE].dropna().unique()
    feature_name = f"on_target_total_hybridization_{cutoff}"
    return _apply_risearch_scoring(
        aso_df=aso_df,
        gene_to_data=gene_to_data,
        target_genes=unique_genes,
        get_gene_fn=lambda row: row[CANONICAL_GENE],
        feature_name=feature_name,
        cutoff=cutoff,
        n_jobs=n_jobs,
        verbose=verbose,
    )


def off_target_specific_seq_pandarallel(aso_df, gene_name, gene_to_data, cutoff, n_jobs=1, verbose=False):
    """Scores each oligo against a single statically provided gene."""
    feature_name = f"off_target_single_{gene_name}_c{cutoff}"
    return _apply_risearch_scoring(
        aso_df=aso_df,
        gene_to_data=gene_to_data,
        target_genes=[gene_name],
        get_gene_fn=lambda row: gene_name,
        feature_name=feature_name,
        cutoff=cutoff,
        n_jobs=n_jobs,
        verbose=verbose,
    )
