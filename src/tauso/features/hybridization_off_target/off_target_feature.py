import logging
import os
import uuid

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


def populate_off_target_specific(ASO_df, gene_to_data, cell_line2data, top_n_list, cutoff_list, method, n_cores=1):
    """
    Enriches ASO_df with off-target scores based on the specific cell line transcriptome.
    Uses CELL_LINE_DEPMAP column to map rows to cell_line2data keys.

    Target FASTA files are built once per (top_n, cell_line) and reused across all cutoffs,
    eliminating redundant disk writes on /dev/shm.
    """
    ASO_df = ASO_df.copy()
    feature_names = []
    grouped = ASO_df.groupby(CELL_LINE_DEPMAP)

    for top_n in top_n_list:
        # Build maps + target files for every cell line once (reused across all cutoffs)
        cell_line_info = {}  # {cell_line: (target_path, seq_map, exp_map)}
        TMP_PATH.mkdir(parents=True, exist_ok=True)

        try:
            for cell_line, _group_df in grouped:
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

            for cutoff in cutoff_list:
                logger.debug("Running Specific Config: TopN=%d, Cutoff=%d", top_n, cutoff)
                spec_col_name = serialize_feature_name(method, top_n, cutoff, is_specific=True)
                calculated_scores = []

                for cell_line, group_df in grouped:
                    if cell_line not in cell_line2data:
                        calculated_scores.append(pd.Series(np.nan, index=group_df.index))
                        continue

                    info = cell_line_info.get(cell_line)
                    if info is None or info[0] is None:
                        calculated_scores.append(pd.Series(0.0, index=group_df.index))
                        continue

                    target_path, spec_seq_map, spec_exp_map = info
                    group_scores = compute_group_batch(
                        group_df,
                        spec_seq_map,
                        spec_exp_map,
                        cutoff,
                        method,
                        prebuilt_target_path=target_path,
                    )
                    calculated_scores.append(group_scores)

                if calculated_scores:
                    full_series = pd.concat(calculated_scores)
                    ASO_df[spec_col_name] = full_series.reindex(ASO_df.index)
                else:
                    ASO_df[spec_col_name] = np.nan

                feature_names.append(spec_col_name)

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
    n_cores=1,
):
    """
    Enriches ASO_df with off-target scores using a single batched RIsearch call per
    (top_n, cutoff) combination instead of one call per row.
    """
    ASO_df = ASO_df.copy()
    feature_names = []

    if "general" not in cell_line2data:
        raise ValueError("Key 'general' not found in cell_line2data dictionary.")
    general_df_all = cell_line2data["general"]

    for top_n in top_n_list:
        general_df = general_df_all.head(top_n)

        # Build seq/exp maps once per top_n
        general_seq_map = {}
        general_exp_map = {}

        norm_col = next((c for c in general_df.columns if "expression_norm" in c), "expression_norm")
        for _, row in general_df.iterrows():
            gene = row["Gene"].split()[0]
            if gene not in gene_to_data:
                raise KeyError(f"CRITICAL: Gene '{gene}' not found in gene_to_data mapping.")
            general_seq_map[gene] = gene_to_data[gene].full_mrna
            general_exp_map[gene] = (
                row.get("expression_TPM", 0),
                row.get(norm_col, row.get("expression_norm", 0)),
            )

        # Build target FASTA once per top_n; reuse it across all cutoffs
        TMP_PATH.mkdir(parents=True, exist_ok=True)
        general_target_path = dump_target_file(f"target-general-{uuid.uuid4().hex}.fa", general_seq_map)
        try:
            for cutoff in cutoff_list:
                logger.debug("Running config: TopN=%d, Cutoff=%d", top_n, cutoff)
                gen_col_name = serialize_feature_name(method, top_n, cutoff, is_specific=False)

                scores = compute_group_batch(
                    ASO_df,
                    general_seq_map,
                    general_exp_map,
                    cutoff,
                    method,
                    prebuilt_target_path=general_target_path,
                )
                ASO_df[gen_col_name] = scores
                feature_names.append(gen_col_name)
        finally:
            if os.path.exists(general_target_path):
                os.remove(general_target_path)

    return ASO_df, feature_names


def populate_off_target_specific_per_rank(
    ASO_df,
    gene_to_data,
    cell_line2data,
    max_rank,
    cutoff_list,
    method,
    n_cores=1,
    verbose=False,
):
    """
    Enriches ASO_df with rank-specific off-target details relative to the specific cell line.

    Instead of summing the Top N genes, this creates columns for the 1st, 2nd, ... Nth
    highest expressed gene individually.

    Returns 3 columns per Rank/Cutoff combo:
      1. Score
      2. Gene Name (String)
      3. Expression Value (Float)
    """
    ASO_df = ASO_df.copy()
    feature_names = []
    accumulator = {}

    logger.info("Grouping by cell line to calculate specific ranks 1 to %d...", max_rank)
    grouped = ASO_df.groupby(CELL_LINE_DEPMAP)

    for cell_line, group_df in grouped:
        if verbose:
            logger.debug("Analyzing cell line: %s", cell_line)
        # 1. Validation: Check if we have data for this cell line
        if cell_line not in cell_line2data:
            logger.warning("Missing cell line: %s", cell_line)
            # Missing context: fill this group with NaNs later
            continue

        # 2. Get the specific transcriptome (only up to the max rank needed)
        specific_df = cell_line2data[cell_line].head(max_rank)

        # 3. Identify Normalization Column (Robustness)
        norm_col = next(
            (c for c in specific_df.columns if "expression_norm" in c),
            "expression_norm",
        )

        # Iterate 0 to max_rank-1 (i.e., Rank 1 to Rank N)
        for rank_idx in range(max_rank):
            if verbose:
                logger.debug("Rank: %d", rank_idx)
            current_rank = rank_idx + 1  # 1-based index for naming

            # If the cell line has fewer genes than the requested rank, skip or fill 0
            if rank_idx >= len(specific_df):
                # You might want to log this or handle it.
                # For now, we just won't compute, leaving it as default (0/NaN) initialization.
                continue

            # Get the gene data for this specific rank
            gene_row = specific_df.iloc[rank_idx]
            gene_name = gene_row["Gene"].split()[0]
            if verbose:
                logger.debug("Analyzing gene: %s", gene_name)
            # Expression value tuple: (TPM, Norm)
            exp_val_tuple = (
                gene_row.get("expression_TPM", 0),
                gene_row.get(norm_col, gene_row.get("expression_norm", 0)),
            )
            # We usually save the normalized value in the dataframe column
            exp_val_scalar = exp_val_tuple[1]

            # Validation: Do we have the sequence?
            if gene_name not in gene_to_data:
                continue

            # Construct Single-Gene Maps for the worker
            single_seq_map = {gene_name: gene_to_data[gene_name].full_mrna}
            single_exp_map = {gene_name: exp_val_tuple}

            for cutoff in cutoff_list:
                # define column names
                base_col = f"OT_Spec_Rank{current_rank}_c{cutoff}"
                col_score = f"{base_col}_Score"
                col_gene = f"{base_col}_Gene"
                col_exp = f"{base_col}_Exp"

                # Ensure lists exist in accumulator
                if col_score not in accumulator:
                    accumulator[col_score] = []
                if col_gene not in accumulator:
                    accumulator[col_gene] = []
                if col_exp not in accumulator:
                    accumulator[col_exp] = []

                scores = compute_group_batch(group_df, single_seq_map, single_exp_map, cutoff, method)

                # Store Results
                accumulator[col_score].append(scores)

                # Create Series for Metadata (Gene Name and Expression)
                # These are constant for the whole group, but must align with index
                gene_series = pd.Series([gene_name] * len(group_df), index=group_df.index)
                exp_series = pd.Series([exp_val_scalar] * len(group_df), index=group_df.index)

                accumulator[col_gene].append(gene_series)
                accumulator[col_exp].append(exp_series)

                # Track feature names (only add unique ones)
                if col_score not in feature_names:
                    feature_names.extend([col_score, col_gene, col_exp])

    # 4. Reassemble Dataframe
    logger.debug("Reassembling chunks...")
    for col_name, chunks in accumulator.items():
        if chunks:
            full_series = pd.concat(chunks)
            # Reindex ensures that rows not processed (e.g. missing cell line data) get NaN/None
            ASO_df[col_name] = full_series.reindex(ASO_df.index)
        else:
            # If no data was computed (empty inputs?), initialize empty
            ASO_df[col_name] = np.nan

    return ASO_df, feature_names
