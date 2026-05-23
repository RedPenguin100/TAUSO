import logging
import os
import uuid
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

from ...data.consts import CANONICAL_GENE, SEQUENCE
from ...mem_debug import log_mem, log_obj_size, track
from ...util import get_antisense
from ..hybridization.fast_hybridization import (
    TMP_PATH,
    Interaction,
    dump_target_file,
    get_triggers_mfe_scores_batch,
    stream_triggers_mfe_hits,
)
from .off_target_functions import parse_risearch_output

_RT = 0.616


def _sum_exp_energy_by_trigger(result_df: pd.DataFrame) -> dict:
    """Return {trigger_id: sum(exp(-RT * energy))} for every trigger in result_df.

    Groups on a raw numpy array (no DataFrame copy). The previous
    `result_df.assign(_exp=...)` cloned the whole frame; for cutoff=0 result
    frames in the hundreds of MB that doubled the transient heap peak.
    """
    exp_vals = np.exp(-_RT * result_df["energy"].to_numpy())
    trigger_arr = result_df["trigger"].to_numpy()
    return pd.Series(exp_vals).groupby(trigger_arr, sort=False).sum().to_dict()


def _validate_genes_found(target_genes, gene_to_data):
    not_found = [g for g in target_genes if g not in gene_to_data]
    if not_found:
        raise ValueError(f"The following genes are not found in gene_to_data: {not_found}")


def _score_one_gene(gene, row_triggers, target_path, cutoff, chunk_size=100, stream=False):
    """Run RIsearch for one gene in chunks and return {trigger_id: score}. Thread-safe.

    Chunking caps peak memory: with cutoff=0 a single 1000-query call can
    generate hundreds of MB of RIsearch output held in RAM before parsing.
    chunk_size=100 gives 10 calls instead of 1 while still being ~100x
    faster than the old one-subprocess-per-row approach.

    stream=True uses the line-streaming RIsearch parser. For cutoff=0 this
    keeps peak Python heap to O(parse_chunk_rows) instead of O(total_hit_bytes),
    skipping the materialized DataFrame entirely.
    """
    if stream:
        return _score_one_gene_streaming(gene, row_triggers, target_path, cutoff, chunk_size)

    combined: dict = {}
    n_chunks = (len(row_triggers) + chunk_size - 1) // chunk_size
    log_mem(f"_score_one_gene[gene={gene}, n_rows={len(row_triggers)}, n_chunks={n_chunks}] start")
    for i in range(0, len(row_triggers), chunk_size):
        chunk = row_triggers[i : i + chunk_size]
        chunk_idx = i // chunk_size
        with track(f"risearch_call[gene={gene}, chunk={chunk_idx}/{n_chunks}, n_q={len(chunk)}]"):
            result = get_triggers_mfe_scores_batch(
                trigger_id_seq_pairs=[(str(idx), trig) for idx, trig in chunk],
                target_file_path=target_path,
                minimum_score=cutoff,
                parsing_type="2",
                interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
                transpose=True,
                batch_id=f"{os.getpid()}-{uuid.uuid4().hex}",
            )
        log_obj_size(f"risearch_result[gene={gene}, chunk={chunk_idx}]", result)
        if not result.strip():
            continue
        with track(f"parse_risearch_output[gene={gene}, chunk={chunk_idx}]"):
            result_df = parse_risearch_output(result)
        log_obj_size(f"result_df[gene={gene}, chunk={chunk_idx}]", result_df)
        del result
        if result_df.empty or "energy" not in result_df.columns:
            continue
        with track(f"_sum_exp_energy[gene={gene}, chunk={chunk_idx}]"):
            combined.update(_sum_exp_energy_by_trigger(result_df))
        del result_df
    log_obj_size(f"combined[gene={gene}]", combined)
    log_mem(f"_score_one_gene[gene={gene}] end")
    return combined


def _score_one_gene_streaming(gene, row_triggers, target_path, cutoff, chunk_size=100):
    """Streaming variant of _score_one_gene: chunk-parses RIsearch stdout via
    pd.read_csv(chunksize) and aggregates exp(-RT·energy) per chunk.

    Avoids materializing the full RIsearch TSV (the giant string + DataFrame
    the non-streaming path holds). Memory is bounded by the pd.read_csv chunk,
    not by total output size — the win for cutoff=0 where output can be
    hundreds of MB per RIsearch call.

    Aggregation is vectorized per chunk (numpy exp + pandas groupby sum) so
    total wall time is competitive with — and on cutoff=0 measurably faster
    than — the materialized path.
    """
    combined: dict = {}
    n_chunks = (len(row_triggers) + chunk_size - 1) // chunk_size
    log_mem(f"_score_one_gene_streaming[gene={gene}, n_rows={len(row_triggers)}, n_chunks={n_chunks}] start")
    for i in range(0, len(row_triggers), chunk_size):
        chunk = row_triggers[i : i + chunk_size]
        chunk_idx = i // chunk_size
        hit_chunk_count = 0
        with track(f"risearch_stream[gene={gene}, chunk={chunk_idx}/{n_chunks}, n_q={len(chunk)}]"):
            for hit_df in stream_triggers_mfe_hits(
                trigger_id_seq_pairs=[(str(idx), trig) for idx, trig in chunk],
                target_file_path=target_path,
                minimum_score=cutoff,
                parsing_type="2",
                interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
                transpose=True,
                batch_id=f"{os.getpid()}-{uuid.uuid4().hex}",
            ):
                hit_chunk_count += 1
                if hit_df.empty:
                    continue
                exp_vals = np.exp(-_RT * hit_df["energy"].to_numpy())
                partial = (
                    pd.Series(exp_vals)
                    .groupby(hit_df["trigger"].to_numpy(), sort=False)
                    .sum()
                )
                # Sum across chunks for the same trigger (rare but possible if a
                # trigger straddles the chunksize boundary).
                for k, v in partial.items():
                    k_str = str(k)
                    existing = combined.get(k_str)
                    combined[k_str] = float(v) if existing is None else existing + float(v)
        log_mem(f"after risearch_stream[gene={gene}, chunk={chunk_idx}, parse_chunks={hit_chunk_count}]")
    log_obj_size(f"combined[gene={gene}]", combined)
    log_mem(f"_score_one_gene_streaming[gene={gene}] end")
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
    stream=False,
):
    """Core logic for RIsearch hybridization scoring."""
    _validate_genes_found(target_genes, gene_to_data)
    TMP_PATH.mkdir(exist_ok=True)

    log_mem(f"_apply_risearch_scoring[feature={feature_name}, n_aso={len(aso_df)}, "
            f"n_genes={len(target_genes)}, n_jobs={n_jobs}, cutoff={cutoff}, stream={stream}] start")

    gene_to_target_info = {}
    try:
        with track(f"build_target_fastas[n_genes={len(target_genes)}]"):
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
            "RIsearch scoring: %d genes × %d rows → %d tasks across %d workers (stream=%s)",
            len(gene_to_row_triggers),
            sum(len(v) for v in gene_to_row_triggers.values()),
            len(tasks),
            effective_workers,
            stream,
        )

        # Accumulate partial results per gene before building the score series
        gene_scores: dict[str, dict] = defaultdict(dict)

        if effective_workers > 1:
            with track(f"thread_pool_dispatch[n_tasks={len(tasks)}, workers={effective_workers}]"):
                with ThreadPoolExecutor(max_workers=effective_workers) as pool:
                    futures = [
                        (
                            pool.submit(
                                _score_one_gene, gene, sub_triggers, target_path, cutoff, stream=stream
                            ),
                            gene,
                        )
                        for gene, sub_triggers, target_path in tasks
                    ]
                    for fut, gene in futures:
                        gene_scores[gene].update(fut.result())
                        log_mem(f"after_task_result[gene={gene}, total_scored={sum(len(d) for d in gene_scores.values())}]")
        else:
            for gene, sub_triggers, target_path in tasks:
                gene_scores[gene].update(
                    _score_one_gene(gene, sub_triggers, target_path, cutoff, stream=stream)
                )

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


def on_target_total_hybridization(aso_df, gene_to_data, cutoff, n_jobs=1, verbose=False, stream=False):
    """Scores each oligo against its own canonical gene.

    stream=True uses line-streaming RIsearch parsing — recommended for cutoff=0
    where the materialized TSV can be hundreds of MB per worker.
    """
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
        stream=stream,
    )


def off_target_specific_seq_pandarallel(
    aso_df, gene_name, gene_to_data, cutoff, n_jobs=1, verbose=False, stream=False
):
    """Scores each oligo against a single statically provided gene.

    stream=True uses line-streaming RIsearch parsing — recommended for cutoff=0
    where the materialized TSV can be hundreds of MB per worker.
    """
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
        stream=stream,
    )
