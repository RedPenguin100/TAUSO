import logging
import os
import uuid
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

from ...data.consts import CELL_LINE_DEPMAP

logger = logging.getLogger(__name__)

from ..hybridization.fast_hybridization import TMP_PATH, dump_target_file
from .add_off_target_feat import compute_group_batch


def serialize_feature_name(method, top_n, cutoff, is_specific):
    if is_specific:
        return f"off_target_score_specific_{method}_n{top_n}_c{cutoff}"
    return f"off_target_score_general_{method}_n{top_n}_c{cutoff}"


def populate_off_target_specific(ASO_df, gene_to_data, cell_line2data, top_n_list, cutoff_list, method, n_jobs=1):
    """
    Enriches ASO_df with off-target scores based on the specific cell line transcriptome.
    Uses CELL_LINE_DEPMAP column to map rows to cell_line2data keys.

    Target FASTA files are built once per (top_n, cell_line) and reused across all cutoffs.
    When n_jobs > 1, cell lines are scored concurrently within each (top_n, cutoff).
    """
    ASO_df = ASO_df.copy()
    feature_names = []
    groups = list(ASO_df.groupby(CELL_LINE_DEPMAP, observed=True))

    for top_n in top_n_list:
        cell_line_info = {}  # {cell_line: (target_path, seq_map, exp_map)}
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
                    cell_line_info[cell_line] = (None, {}, {})
                    continue

                target_path = dump_target_file(f"target-spec-{cell_line}-{uuid.uuid4().hex}.fa", spec_seq_map)
                cell_line_info[cell_line] = (target_path, spec_seq_map, spec_exp_map)

            def _score_cell_line(cell_line, group_df, cutoff, _info_map):
                if cell_line not in cell_line2data:
                    return pd.Series(np.nan, index=group_df.index)
                info = _info_map.get(cell_line)
                if info is None or info[0] is None:
                    return pd.Series(0.0, index=group_df.index)
                target_path, seq_map, exp_map = info
                return compute_group_batch(
                    group_df, seq_map, exp_map, cutoff, method, prebuilt_target_path=target_path
                )

            for cutoff in cutoff_list:
                logger.debug("Running Specific Config: TopN=%d, Cutoff=%d", top_n, cutoff)
                col = serialize_feature_name(method, top_n, cutoff, is_specific=True)

                if n_jobs > 1 and len(groups) > 1:
                    with ThreadPoolExecutor(max_workers=min(n_jobs, len(groups))) as pool:
                        future_to_idx = {
                            pool.submit(_score_cell_line, cl, gdf, cutoff, cell_line_info): i
                            for i, (cl, gdf) in enumerate(groups)
                        }
                        ordered = [None] * len(groups)
                        for fut, idx in future_to_idx.items():
                            ordered[idx] = fut.result()
                    full_series = pd.concat(ordered)
                else:
                    full_series = pd.concat(
                        [_score_cell_line(cl, gdf, cutoff, cell_line_info) for cl, gdf in groups]
                    )

                ASO_df[col] = full_series.reindex(ASO_df.index)
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
):
    """
    Enriches ASO_df with off-target scores using a single batched RIsearch call per
    (top_n, cutoff) combination.

    Target FASCTAs are built once per top_n, then all (top_n, cutoff) combinations
    are run — concurrently when n_jobs > 1. Results are applied in original order.
    """
    ASO_df = ASO_df.copy()
    feature_names = []

    if "general" not in cell_line2data:
        raise ValueError("Key 'general' not found in cell_line2data dictionary.")
    general_df_all = cell_line2data["general"]

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

        # Phase 2: run all (top_n, cutoff) combinations, optionally in parallel
        combos = [(top_n, cutoff) for top_n in top_n_list for cutoff in cutoff_list]

        def _run_combo(top_n, cutoff):
            seq_map, exp_map, target_path = top_n_data[top_n]
            col = serialize_feature_name(method, top_n, cutoff, is_specific=False)
            logger.debug("Running config: TopN=%d, Cutoff=%d", top_n, cutoff)
            scores = compute_group_batch(
                ASO_df, seq_map, exp_map, cutoff, method, prebuilt_target_path=target_path
            )
            return col, scores

        if n_jobs > 1 and len(combos) > 1:
            results = {}
            with ThreadPoolExecutor(max_workers=min(n_jobs, len(combos))) as pool:
                future_to_combo = {pool.submit(_run_combo, tn, c): (tn, c) for tn, c in combos}
                for fut, combo in future_to_combo.items():
                    results[combo] = fut.result()
            for combo in combos:
                col, scores = results[combo]
                ASO_df[col] = scores
                feature_names.append(col)
        else:
            for top_n, cutoff in combos:
                col, scores = _run_combo(top_n, cutoff)
                ASO_df[col] = scores
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

                scores = compute_group_batch(group_df, single_seq_map, single_exp_map, cutoff, method)
                partial.setdefault(col_score, []).append(scores)
                partial.setdefault(col_gene, []).append(
                    pd.Series([gene_name] * len(group_df), index=group_df.index)
                )
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
