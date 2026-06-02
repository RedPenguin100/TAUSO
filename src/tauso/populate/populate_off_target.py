import logging
import os
import uuid

import numpy as np
import pandas as pd

from ..data.consts import CELL_LINE_DEPMAP
from ..features.hybridization.fast_hybridization import TMP_PATH, dump_target_file
from ..features.hybridization.off_target.add_off_target_feat import compute_group_batch_multi_cutoff_multi_topn
from ..features.hybridization.off_target.parallel import run_tasks_parallel

logger = logging.getLogger(__name__)


def serialize_feature_name(method, top_n, cutoff, is_specific):
    if is_specific:
        return f"off_target_score_specific_{method}_n{top_n}_c{cutoff}"
    return f"off_target_score_general_{method}_n{top_n}_c{cutoff}"


def _chunk_df(df, chunk_size):
    """Split df into a list of sub-DataFrames of at most chunk_size rows."""
    return [df.iloc[i : i + chunk_size] for i in range(0, len(df), chunk_size)]


def _score_chunk_multi_cutoff_multi_topn(chunk_df, top_n_to_data, cutoffs, method, target_path):
    """One RIsearch pass for a chunk against a max(top_n) target, deriving per-top_n
    results by filtering to each top_n's gene_set. Returns {(top_n, cutoff): Series}."""
    return compute_group_batch_multi_cutoff_multi_topn(
        chunk_df, top_n_to_data, cutoffs, method, prebuilt_target_path=target_path
    )


def _build_top_n_to_data(expression_df, top_n_list, gene_to_data, *, require_seq=True):
    """For each top_n, return (exp_map, gene_set) for head(top_n) of expression_df.

    Also returns the union seq_map keyed by gene → mRNA sequence, for genes in the
    largest gene_set that have entries in gene_to_data. If ``require_seq`` is True
    (general path), a missing gene_to_data entry is fatal; if False (specific path),
    such genes are silently skipped (matching the previous behaviour).
    """
    top_n_to_data: dict = {}
    seq_map: dict = {}
    max_top_n = max(top_n_list) if top_n_list else 0
    union_df = expression_df.head(max_top_n)
    norm_col = next((c for c in union_df.columns if "expression_norm" in c), "expression_norm")
    for top_n in sorted(set(top_n_list)):
        sub = expression_df.head(top_n)
        exp_map: dict = {}
        gene_set: set = set()
        for _, row in sub.iterrows():
            gene = row["Gene"].split()[0]
            if gene not in gene_to_data:
                if require_seq:
                    raise KeyError(f"CRITICAL: Gene '{gene}' not found in gene_to_data mapping.")
                continue
            exp_map[gene] = (
                row.get("expression_TPM", 0),
                row.get(norm_col, row.get("expression_norm", 0)),
            )
            gene_set.add(gene)
        top_n_to_data[top_n] = (exp_map, gene_set)
    # Build seq_map for the union (all genes in any gene_set).
    union_genes: set = set()
    for _, gene_set in top_n_to_data.values():
        union_genes |= gene_set
    for gene in union_genes:
        seq_map[gene] = gene_to_data[gene].full_mrna
    return top_n_to_data, seq_map


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

    Top_n-collapse and cutoff-collapse applied per cell line: one RIsearch pass per
    (cell_line, ASO chunk) at the loosest cutoff against the head(max(top_n_list))
    target, every (top_n, cutoff) feature derived from it (smaller top_n by gene-subset
    filter, cutoffs by score-filter on the streaming pyarrow output). ASO chunks across
    all cell lines share one thread pool, so the full run already saturates the workers
    (no gene-sharding).

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

    TMP_PATH.mkdir(parents=True, exist_ok=True)
    # cell_info[cell_line] = {is_known, top_n_to_data, target_path, chunks}
    cell_info: dict = {}
    created_paths: list = []
    try:
        for cell_line, group_df in groups:
            chunks = _chunk_df(group_df, chunk_size)
            if cell_line not in cell_line2data:
                cell_info[cell_line] = {"is_known": False, "target_path": None, "chunks": chunks}
                continue

            top_n_to_data, seq_map = _build_top_n_to_data(
                cell_line2data[cell_line], top_n_list, gene_to_data, require_seq=False
            )

            target_path = None
            if seq_map:
                target_path = dump_target_file(f"target-spec-{cell_line}-{uuid.uuid4().hex}.fa", seq_map)
                created_paths.append(target_path)
            cell_info[cell_line] = {
                "is_known": True,
                "top_n_to_data": top_n_to_data,
                "target_path": target_path,
                "chunks": chunks,
            }

        tasks = [
            ((cell_line, chunk_idx), chunk_df, info["top_n_to_data"], cutoff_list, method, info["target_path"])
            for cell_line, info in cell_info.items()
            if info["is_known"] and info["target_path"] is not None
            for chunk_idx, chunk_df in enumerate(info["chunks"])
        ]
        results = run_tasks_parallel(tasks, _score_chunk_multi_cutoff_multi_topn, n_jobs)

        for top_n in top_n_list:
            for cutoff in cutoff_list:
                col = serialize_feature_name(method, top_n, cutoff, is_specific=True)
                series_list = []
                for cell_line, info in cell_info.items():
                    if not info["is_known"]:  # unknown cell line -> NaN
                        series_list += [pd.Series(np.nan, index=c.index) for c in info["chunks"]]
                    elif info["target_path"] is None:  # known, no usable target -> 0.0
                        series_list += [pd.Series(0.0, index=c.index) for c in info["chunks"]]
                    else:
                        series_list += [results[(cell_line, i)][(top_n, cutoff)] for i in range(len(info["chunks"]))]
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

    Top_n-collapse and cutoff-collapse: one RIsearch pass per ASO chunk at the loosest
    cutoff against the head(max(top_n_list)) target, every (top_n, cutoff) feature
    derived from it (smaller top_n by gene-subset filter, cutoffs by score-filter on
    the streaming pyarrow output). Each chunk is at most chunk_size ASOs to bound peak
    RIsearch memory; results are assembled in original (top_n, cutoff) order. `stream`
    is accepted for backward compatibility; the multi-cutoff path always streams.
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

    TMP_PATH.mkdir(parents=True, exist_ok=True)
    top_n_to_data, seq_map = _build_top_n_to_data(general_df_all, top_n_list, gene_to_data, require_seq=True)
    target_path = dump_target_file(f"target-general-{uuid.uuid4().hex}.fa", seq_map)

    try:
        tasks = [
            ((chunk_idx,), chunk_df, top_n_to_data, cutoff_list, method, target_path)
            for chunk_idx, chunk_df in enumerate(aso_chunks)
        ]
        # results[(chunk_idx,)] = {(top_n, cutoff): Series}
        results = run_tasks_parallel(tasks, _score_chunk_multi_cutoff_multi_topn, n_jobs)

        for top_n in top_n_list:
            for cutoff in cutoff_list:
                col = serialize_feature_name(method, top_n, cutoff, is_specific=False)
                full_series = pd.concat([results[(i,)][(top_n, cutoff)] for i in range(len(aso_chunks))])
                ASO_df[col] = full_series.reindex(ASO_df.index)
                feature_names.append(col)

    finally:
        if os.path.exists(target_path):
            os.remove(target_path)

    return ASO_df, feature_names
