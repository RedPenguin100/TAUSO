import logging
import os
import uuid
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

from ..data.consts import CELL_LINE_DEPMAP
from ..features.hybridization.fast_hybridization import TMP_PATH, dump_target_file
from ..features.hybridization.off_target.add_off_target_feat import compute_group_batch_multi_cutoff

logger = logging.getLogger(__name__)


def serialize_feature_name(method, top_n, cutoff, is_specific):
    if is_specific:
        return f"off_target_score_specific_{method}_n{top_n}_c{cutoff}"
    return f"off_target_score_general_{method}_n{top_n}_c{cutoff}"


def _chunk_df(df, chunk_size):
    """Split df into a list of sub-DataFrames of at most chunk_size rows."""
    return [df.iloc[i : i + chunk_size] for i in range(0, len(df), chunk_size)]


def _score_chunk_multi_cutoff(chunk_df, exp_map, cutoffs, method, target_path):
    """One RIsearch pass for a chunk against a prebuilt target, scoring every cutoff
    from that single pass. Returns {cutoff: Series}."""
    return compute_group_batch_multi_cutoff(chunk_df, exp_map, cutoffs, method, prebuilt_target_path=target_path)


def _run_tasks_parallel(tasks, fn, n_jobs):
    """
    Run a list of (key, *args) tasks using fn(*args).
    Returns {key: result} in submission order.
    When n_jobs == 1, runs serially.
    """
    results = {}
    if n_jobs > 1 and len(tasks) > 1:
        with ThreadPoolExecutor(max_workers=min(n_jobs, len(tasks))) as pool:
            future_to_key = {pool.submit(fn, *args): key for key, *args in tasks}
            for fut, key in future_to_key.items():
                results[key] = fut.result()
    else:
        for key, *args in tasks:
            results[key] = fn(*args)
    return results


def populate_off_target_specific(
    ASO_df,
    gene_to_data,
    cell_line2data,
    top_n_list,
    cutoff_list,
    method,
    n_jobs=1,
    chunk_size=250,
    stream=True,
):
    """
    Enriches ASO_df with off-target scores based on the specific cell line transcriptome.

    Cutoff-collapse applied per cell line: one RIsearch pass per (cell_line, ASO chunk)
    at the loosest cutoff, every cutoff derived from it. ASO chunks across all cell lines
    share one thread pool, so the full run already saturates the workers (no gene-sharding).

    Fallbacks: a cell line absent from cell_line2data scores NaN; a known cell line with no
    usable target genes scores 0.0. `stream` is accepted for backward compatibility; the
    multi-cutoff path always streams.
    """
    ASO_df = ASO_df.copy()
    feature_names = []
    groups = list(ASO_df.groupby(CELL_LINE_DEPMAP, observed=True))

    logger.info(
        "populate_off_target_specific: top_n=%s cutoffs=%s n_aso=%d n_cell_lines=%d n_jobs=%d",
        top_n_list,
        cutoff_list,
        len(ASO_df),
        len(groups),
        n_jobs,
    )

    for top_n in top_n_list:
        TMP_PATH.mkdir(parents=True, exist_ok=True)
        # cell_info[cell_line] = {is_known, exp_map, target_path, chunks}
        cell_info: dict = {}
        created_paths = []
        try:
            for cell_line, group_df in groups:
                chunks = _chunk_df(group_df, chunk_size)
                if cell_line not in cell_line2data:
                    cell_info[cell_line] = {"is_known": False, "target_path": None, "chunks": chunks}
                    continue

                specific_df = cell_line2data[cell_line].head(top_n)
                norm_col = next((c for c in specific_df.columns if "expression_norm" in c), "expression_norm")
                spec_seq_map: dict = {}
                spec_exp_map: dict = {}
                for _, row in specific_df.iterrows():
                    gene = row["Gene"].split()[0]
                    if gene in gene_to_data:
                        spec_seq_map[gene] = gene_to_data[gene].full_mrna
                        spec_exp_map[gene] = (
                            row.get("expression_TPM", 0),
                            row.get(norm_col, row.get("expression_norm", 0)),
                        )

                target_path = None
                if spec_seq_map:
                    target_path = dump_target_file(f"target-spec-{cell_line}-{uuid.uuid4().hex}.fa", spec_seq_map)
                    created_paths.append(target_path)
                cell_info[cell_line] = {
                    "is_known": True,
                    "exp_map": spec_exp_map,
                    "target_path": target_path,
                    "chunks": chunks,
                }

            # One RIsearch pass per (cell_line, chunk); all cutoffs derived from it.
            tasks = [
                ((cell_line, chunk_idx), chunk_df, info["exp_map"], cutoff_list, method, info["target_path"])
                for cell_line, info in cell_info.items()
                if info["is_known"] and info["target_path"] is not None
                for chunk_idx, chunk_df in enumerate(info["chunks"])
            ]
            results = _run_tasks_parallel(tasks, _score_chunk_multi_cutoff, n_jobs)

            for cutoff in cutoff_list:
                col = serialize_feature_name(method, top_n, cutoff, is_specific=True)
                series_list = []
                for cell_line, info in cell_info.items():
                    if not info["is_known"]:  # unknown cell line -> NaN
                        series_list += [pd.Series(np.nan, index=c.index) for c in info["chunks"]]
                    elif info["target_path"] is None:  # known, no usable target -> 0.0
                        series_list += [pd.Series(0.0, index=c.index) for c in info["chunks"]]
                    else:
                        series_list += [results[(cell_line, i)][cutoff] for i in range(len(info["chunks"]))]
                ASO_df[col] = pd.concat(series_list).reindex(ASO_df.index)
                feature_names.append(col)

        finally:
            for path in created_paths:
                if os.path.exists(path):
                    os.remove(path)

    return ASO_df, feature_names


def populate_off_target_general(
    ASO_df,
    gene_to_data,
    cell_line2data,
    top_n_list,
    cutoff_list,
    method,
    n_jobs=1,
    chunk_size=250,
    stream=True,
):
    """
    Enriches ASO_df with off-target scores using batched RIsearch calls.

    One RIsearch pass per (top_n, ASO chunk) runs at the loosest cutoff, and every
    cutoff is derived from that single streaming pass (cutoff-collapse). The full
    173k-ASO run produces ~700 chunk tasks, which already saturate the workers, so
    there is no gene-sharding here. Each chunk is at most chunk_size ASOs to bound
    peak RIsearch memory; results are assembled in original (top_n, cutoff) order.
    `stream` is accepted for backward compatibility; the multi-cutoff path always
    streams.
    """
    ASO_df = ASO_df.copy()
    feature_names = []

    if "general" not in cell_line2data:
        raise ValueError("Key 'general' not found in cell_line2data dictionary.")
    general_df_all = cell_line2data["general"]

    aso_chunks = _chunk_df(ASO_df, chunk_size)

    logger.info(
        "populate_off_target_general: top_n=%s cutoffs=%s n_aso=%d n_chunks=%d n_jobs=%d",
        top_n_list,
        cutoff_list,
        len(ASO_df),
        len(aso_chunks),
        n_jobs,
    )

    # Phase 1: build one target FASTA per top_n.
    TMP_PATH.mkdir(parents=True, exist_ok=True)
    top_n_data = {}  # {top_n: (exp_map, target_path)}
    try:
        for top_n in top_n_list:
            general_df = general_df_all.head(top_n)
            seq_map: dict = {}
            exp_map: dict = {}
            norm_col = next((c for c in general_df.columns if "expression_norm" in c), "expression_norm")
            for _, row in general_df.iterrows():
                gene = row["Gene"].split()[0]
                if gene not in gene_to_data:
                    raise KeyError(f"CRITICAL: Gene '{gene}' not found in gene_to_data mapping.")
                seq_map[gene] = gene_to_data[gene].full_mrna
                exp_map[gene] = (
                    row.get("expression_TPM", 0),
                    row.get(norm_col, row.get("expression_norm", 0)),
                )
            target_path = dump_target_file(f"target-general-{uuid.uuid4().hex}.fa", seq_map)
            top_n_data[top_n] = (exp_map, target_path)

        # Phase 2: one RIsearch pass per (top_n, chunk); all cutoffs derived from it.
        tasks = [
            ((top_n, chunk_idx), chunk_df, top_n_data[top_n][0], cutoff_list, method, top_n_data[top_n][1])
            for top_n in top_n_list
            for chunk_idx, chunk_df in enumerate(aso_chunks)
        ]
        # results[(top_n, chunk_idx)] = {cutoff: Series}
        results = _run_tasks_parallel(tasks, _score_chunk_multi_cutoff, n_jobs)

        for top_n in top_n_list:
            for cutoff in cutoff_list:
                col = serialize_feature_name(method, top_n, cutoff, is_specific=False)
                full_series = pd.concat([results[(top_n, i)][cutoff] for i in range(len(aso_chunks))])
                ASO_df[col] = full_series.reindex(ASO_df.index)
                feature_names.append(col)

    finally:
        for _exp_map, target_path in top_n_data.values():
            if os.path.exists(target_path):
                os.remove(target_path)

    return ASO_df, feature_names
