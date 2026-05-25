import logging
import os
import uuid
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

import numpy as np
import pandas as pd

from ..data.consts import CELL_LINE_DEPMAP
from ..features.hybridization.fast_hybridization import TMP_PATH, dump_target_file
from ..features.hybridization.off_target.add_off_target_feat import compute_group_batch

logger = logging.getLogger(__name__)


def _select_pool(cutoff_list):
    """Pool for a batch of scoring tasks. cutoff=0 is parse-bound and the parse
    step is GIL-bound, so a process pool scales where threads plateau; higher
    cutoffs emit little output so threads (lower overhead) are best.
    TAUSO_OFFTARGET_POOL=thread|process overrides.
    """
    override = os.environ.get("TAUSO_OFFTARGET_POOL")
    if override == "process":
        return ProcessPoolExecutor
    if override == "thread":
        return ThreadPoolExecutor
    return ProcessPoolExecutor if 0 in cutoff_list else ThreadPoolExecutor


def serialize_feature_name(method, top_n, cutoff, is_specific):
    if is_specific:
        return f"off_target_score_specific_{method}_n{top_n}_c{cutoff}"
    return f"off_target_score_general_{method}_n{top_n}_c{cutoff}"


def _chunk_df(df, chunk_size):
    """Split df into a list of sub-DataFrames of at most chunk_size rows."""
    return [df.iloc[i : i + chunk_size] for i in range(0, len(df), chunk_size)]


def _score_chunk(chunk_df, exp_map, cutoff, method, target_path, stream=True):
    """Single atomic RIsearch scoring unit. Thread/process-safe: uses a unique query
    FASTA per call and reads the prebuilt target from `target_path` (the seq_map is
    only needed to build that file, which the caller already did)."""
    return compute_group_batch(chunk_df, None, exp_map, cutoff, method, prebuilt_target_path=target_path, stream=stream)


def _score_chunk_or_fallback(chunk_df, info, is_known_cell_line, cutoff, method, stream=True):
    """Like _score_chunk but handles missing/empty cell-line info."""
    if not is_known_cell_line:
        return pd.Series(np.nan, index=chunk_df.index)
    if info is None or info[0] is None:
        return pd.Series(0.0, index=chunk_df.index)
    target_path, exp_map = info
    return compute_group_batch(chunk_df, None, exp_map, cutoff, method, prebuilt_target_path=target_path, stream=stream)


def _run_tasks_parallel(tasks, fn, n_jobs, pool_cls=ThreadPoolExecutor):
    """
    Run a list of (key, *args) tasks using fn(*args).
    Returns {key: result} in submission order.
    When n_jobs == 1, runs serially.
    """
    results = {}
    if n_jobs > 1 and len(tasks) > 1:
        with pool_cls(max_workers=min(n_jobs, len(tasks))) as pool:
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

    Target FASCTAs are built once per (top_n, cell_line) and reused across all cutoffs.
    Parallelism: all (cutoff, cell_line, aso_chunk) combinations run concurrently up to
    n_jobs threads, with each chunk capped at chunk_size ASOs to bound peak RIsearch memory.

    stream=True (default) uses the streaming RIsearch parser inside compute_group_batch
    so per-task peak memory stays bounded regardless of cutoff.
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
        cell_line_info = {}
        TMP_PATH.mkdir(parents=True, exist_ok=True)

        try:
            for cell_line, _group_df in groups:
                if cell_line not in cell_line2data:
                    continue

                specific_df = cell_line2data[cell_line].head(top_n)
                norm_col = next(
                    (c for c in specific_df.columns if "expression_norm" in c),
                    "expression_norm",
                )

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

                if not spec_seq_map:
                    cell_line_info[cell_line] = (None, {})
                    continue

                target_path = dump_target_file(f"target-spec-{cell_line}-{uuid.uuid4().hex}.fa", spec_seq_map)
                cell_line_info[cell_line] = (target_path, spec_exp_map)

            # Pre-chunk every cell line's group so chunks are shared across cutoffs
            cell_line_chunks = {cl: _chunk_df(gdf, chunk_size) for cl, gdf in groups}

            # Tasks: (key, chunk_df, info, is_known, cutoff, method, stream)
            # key = (cutoff, cell_line, chunk_idx) for ordered reassembly
            tasks = [
                (
                    (cutoff, cell_line, chunk_idx),
                    chunk_df,
                    cell_line_info.get(cell_line),
                    cell_line in cell_line2data,
                    cutoff,
                    method,
                    stream,
                )
                for cutoff in cutoff_list
                for cell_line, _ in groups
                for chunk_idx, chunk_df in enumerate(cell_line_chunks[cell_line])
            ]

            n_threads = min(n_jobs, len(tasks))
            logger.debug(
                "populate_off_target_specific(top_n=%d): %d tasks "
                "(%d cutoffs × %d cell_lines × chunks≤%d), using %d thread(s)",
                top_n,
                len(tasks),
                len(cutoff_list),
                len(groups),
                chunk_size,
                n_threads,
            )
            results = _run_tasks_parallel(tasks, _score_chunk_or_fallback, n_jobs, pool_cls=_select_pool(cutoff_list))

            for cutoff in cutoff_list:
                col = serialize_feature_name(method, top_n, cutoff, is_specific=True)

                series_list = []
                for cell_line, _ in groups:
                    n_chunks = len(cell_line_chunks[cell_line])
                    cell_series = pd.concat([results[(cutoff, cell_line, i)] for i in range(n_chunks)])
                    series_list.append(cell_series)

                ASO_df[col] = pd.concat(series_list).reindex(ASO_df.index)
                feature_names.append(col)

        finally:
            for info in cell_line_info.values():
                path = info[0]
                if path is not None and os.path.exists(path):
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
    Enriches ASO_df with off-target scores using a batched RIsearch call per chunk.

    Target FASCTAs are built once per top_n, then all (top_n, cutoff, aso_chunk)
    combinations are run — concurrently when n_jobs > 1. Each chunk is at most
    chunk_size ASOs to bound peak RIsearch memory. Results are assembled in
    original (top_n, cutoff) order.

    stream=True (default) uses the streaming RIsearch parser inside
    compute_group_batch so per-task peak memory stays bounded regardless of cutoff.
    """
    ASO_df = ASO_df.copy()
    feature_names = []

    if "general" not in cell_line2data:
        raise ValueError("Key 'general' not found in cell_line2data dictionary.")
    general_df_all = cell_line2data["general"]

    logger.info(
        "populate_off_target_general: top_n=%s cutoffs=%s n_aso=%d n_jobs=%d",
        top_n_list,
        cutoff_list,
        len(ASO_df),
        n_jobs,
    )

    aso_chunks = _chunk_df(ASO_df, chunk_size)

    # Phase 1: build all target FASCTAs (one per top_n) — serial, disk I/O
    TMP_PATH.mkdir(parents=True, exist_ok=True)
    top_n_data = {}  # {top_n: (seq_map, exp_map, target_path)}
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
            top_n_data[top_n] = (seq_map, exp_map, target_path)

        # Phase 2: tasks = all (top_n, cutoff, chunk_idx) combinations
        tasks = [
            (
                (top_n, cutoff, chunk_idx),
                chunk_df,
                top_n_data[top_n][1],  # exp_map
                cutoff,
                method,
                top_n_data[top_n][2],  # prebuilt target_path
                stream,
            )
            for top_n in top_n_list
            for cutoff in cutoff_list
            for chunk_idx, chunk_df in enumerate(aso_chunks)
        ]

        if n_jobs > 1:
            logger.debug("Dispatching %d tasks across %d threads", len(tasks), min(n_jobs, len(tasks)))

        results = _run_tasks_parallel(tasks, _score_chunk, n_jobs, pool_cls=_select_pool(cutoff_list))

        for top_n in top_n_list:
            for cutoff in cutoff_list:
                col = serialize_feature_name(method, top_n, cutoff, is_specific=False)
                full_series = pd.concat([results[(top_n, cutoff, i)] for i in range(len(aso_chunks))])
                ASO_df[col] = full_series.reindex(ASO_df.index)
                feature_names.append(col)

    finally:
        for _, _, target_path in top_n_data.values():
            if os.path.exists(target_path):
                os.remove(target_path)

    return ASO_df, feature_names


def populate_off_target_specific_per_rank(
    ASO_df,
    gene_to_data,
    cell_line2data,
    max_rank,
    cutoff_list,
    method,
    n_jobs=1,
    stream=True,
):
    """
    Enriches ASO_df with rank-specific off-target details relative to the specific cell line.

    Instead of summing the Top N genes, this creates columns for the 1st, 2nd, ... Nth
    highest expressed gene individually.

    Returns 3 columns per Rank/Cutoff combo:
      1. Score
      2. Gene Name (String)
      3. Expression Value (Float)

    When n_jobs > 1, cell lines are processed concurrently.

    stream=True (default) uses the streaming RIsearch parser inside
    compute_group_batch so per-task peak memory stays bounded regardless of cutoff.
    """
    ASO_df = ASO_df.copy()
    feature_names = []
    accumulator: dict = {}

    logger.info("Grouping by cell line to calculate specific ranks 1 to %d...", max_rank)
    groups = list(ASO_df.groupby(CELL_LINE_DEPMAP, observed=True))

    def _process_cell_line(cell_line, group_df):
        """Returns partial {col_name: series} for one cell line."""
        partial: dict = {}
        if cell_line not in cell_line2data:
            logger.warning("Warning! missing cell line: %s", cell_line)
            return partial

        specific_df = cell_line2data[cell_line].head(max_rank)
        norm_col = next(
            (c for c in specific_df.columns if "expression_norm" in c),
            "expression_norm",
        )

        for rank_idx in range(max_rank):
            logger.debug("Rank: %d", rank_idx)
            if rank_idx >= len(specific_df):
                continue

            gene_row = specific_df.iloc[rank_idx]
            gene_name = gene_row["Gene"].split()[0]
            logger.debug("Analyzing gene: %s", gene_name)

            if gene_name not in gene_to_data:
                continue

            exp_val_tuple = (
                gene_row.get("expression_TPM", 0),
                gene_row.get(norm_col, gene_row.get("expression_norm", 0)),
            )
            exp_val_scalar = exp_val_tuple[1]
            single_seq_map = {gene_name: gene_to_data[gene_name].full_mrna}
            single_exp_map = {gene_name: exp_val_tuple}
            current_rank = rank_idx + 1

            for cutoff in cutoff_list:
                base_col = f"OT_Spec_Rank{current_rank}_c{cutoff}"
                col_score = f"{base_col}_Score"
                col_gene = f"{base_col}_Gene"
                col_exp = f"{base_col}_Exp"

                scores = compute_group_batch(group_df, single_seq_map, single_exp_map, cutoff, method, stream=stream)
                partial.setdefault(col_score, []).append(scores)
                partial.setdefault(col_gene, []).append(pd.Series([gene_name] * len(group_df), index=group_df.index))
                partial.setdefault(col_exp, []).append(
                    pd.Series([exp_val_scalar] * len(group_df), index=group_df.index)
                )

        return partial

    if n_jobs > 1 and len(groups) > 1:
        with ThreadPoolExecutor(max_workers=min(n_jobs, len(groups))) as pool:
            futures = [pool.submit(_process_cell_line, cl, gdf) for cl, gdf in groups]
            for fut in futures:
                for col, series_list in fut.result().items():
                    accumulator.setdefault(col, []).extend(series_list)
    else:
        for cell_line, group_df in groups:
            for col, series_list in _process_cell_line(cell_line, group_df).items():
                accumulator.setdefault(col, []).extend(series_list)

    # Track feature names in rank/cutoff order (not cell_line order)
    for rank_idx in range(max_rank):
        current_rank = rank_idx + 1
        for cutoff in cutoff_list:
            base_col = f"OT_Spec_Rank{current_rank}_c{cutoff}"
            for suffix in ("_Score", "_Gene", "_Exp"):
                col = base_col + suffix
                if col in accumulator:
                    feature_names.append(col)

    logger.debug("Reassembling chunks...")
    for col_name, chunks in accumulator.items():
        if chunks:
            ASO_df[col_name] = pd.concat(chunks).reindex(ASO_df.index)
        else:
            ASO_df[col_name] = np.nan

    return ASO_df, feature_names
