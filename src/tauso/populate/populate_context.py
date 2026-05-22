import logging
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from ..data.consts import CANONICAL_GENE, CELL_LINE_DEPMAP
from ..features.context.ribo_seq import (
    feature_names,
    get_ribo_40s_human_data,
    process_gene_group,
)

logger = logging.getLogger(__name__)


def _run_gene_groups_parallel(bw_path, gene_groups, flanks, how, n_jobs):
    """
    Dispatch process_gene_group over all (gene, gene_rows) pairs.

    Each worker opens its own BigWig handle so threads don't share state.
    Returns (results, skipped) where results is {df_index: feat_dict} and
    skipped is {contig_str: set_of_gene_names}.
    """
    results = {}
    skipped = {}

    def _collect(gene, gene_rows):
        try:
            return gene, process_gene_group(bw_path, gene_rows, flanks, how), None
        except KeyError as exc:
            return gene, {}, str(exc).strip("'\"")

    if n_jobs > 1 and len(gene_groups) > 1:
        with ThreadPoolExecutor(max_workers=min(n_jobs, len(gene_groups))) as pool:
            futures = {pool.submit(_collect, gene, rows): gene for gene, rows in gene_groups}
            for fut in futures:
                _, group_results, contig_err = fut.result()
                gene = futures[fut]
                if contig_err is not None:
                    skipped.setdefault(contig_err, set()).add(gene)
                else:
                    results.update(group_results)
    else:
        for gene, gene_rows in gene_groups:
            _, group_results, contig_err = _collect(gene, gene_rows)
            if contig_err is not None:
                skipped.setdefault(contig_err, set()).add(gene)
            else:
                results.update(group_results)

    return results, skipped


def populate_ribo_seq(organism, aso_df, flanks=(0, 10, 20, 50, 100, 125, 150), how="mean", n_jobs=1):
    if organism != "human":
        raise ValueError("Unsupported organism for ribo_seq feature")
    if how not in {"sum", "mean", "max", "nz_mean"}:
        raise ValueError(how)

    bw_path = str(get_ribo_40s_human_data())
    feat_cols = feature_names(flanks, how)

    for col in feat_cols:
        aso_df[col] = np.nan

    valid_mask = aso_df["chrom"].notna() & aso_df["target_start"].notna()
    valid_df = aso_df[valid_mask]
    n_valid = len(valid_df)
    n_genes = valid_df[CANONICAL_GENE].nunique()

    logger.info(
        "populate_ribo_seq: %d valid rows across %d genes × %d features (n_jobs=%d)",
        n_valid, n_genes, len(feat_cols), n_jobs,
    )

    gene_groups = list(valid_df.groupby(CANONICAL_GENE))

    t0 = time.perf_counter()
    results, skipped = _run_gene_groups_parallel(bw_path, gene_groups, flanks, how, n_jobs)
    elapsed = time.perf_counter() - t0

    logger.info(
        "populate_ribo_seq: done in %.2fs (%.0f rows/s)",
        elapsed,
        n_valid / elapsed if elapsed > 0 else 0,
    )

    if results:
        result_df = pd.DataFrame.from_dict(results, orient="index")
        aso_df.update(result_df)

    for contig, genes in skipped.items():
        logger.warning(
            "Skipped ribo-seq for contig '%s' (genes: %s). Rows will have NaN features.",
            contig,
            ", ".join(sorted(genes)),
        )

    return aso_df, feat_cols


def populate_mrna_expression(
    df: pd.DataFrame, expression_dict: Dict[str, pd.DataFrame]
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Enriches the main dataframe with mRNA expression data for the target gene
    and the specific RNASEH1 gene.
    """

    # --- MINIMAL FIX 1: Warn and drop to prevent column duplication ---
    existing_cols = [c for c in ["target_expression", "rnase_expression"] if c in df.columns]
    if existing_cols:
        logger.warning("Dropping existing columns to prevent alignment breaks: %s", existing_cols)
        df = df.drop(columns=existing_cols)
    # ------------------------------------------------------------------

    # 1. Flatten the dictionary into one Master DataFrame
    # ---------------------------------------------------------
    dfs_to_concat = []

    # We iterate through the dictionary to add the cell line key as a column
    for depmap_id, t_df in expression_dict.items():
        # Defensive copy to avoid SettingWithCopy warnings on the source DFs
        temp = t_df[["Gene", "expression_norm"]].copy()
        temp[CELL_LINE_DEPMAP] = depmap_id
        dfs_to_concat.append(temp)

    expression_master = pd.concat(dfs_to_concat, ignore_index=True)

    # --- MINIMAL FIX 2: Prevent row multiplication during merge ---
    expression_master = expression_master.drop_duplicates(subset=[CELL_LINE_DEPMAP, "Gene"])
    # --------------------------------------------------------------

    # 2. Get 'target_expression' (General Merge)
    # ---------------------------------------------------------
    enhanced_df = df.merge(
        expression_master,
        left_on=[CELL_LINE_DEPMAP, CANONICAL_GENE],
        right_on=[CELL_LINE_DEPMAP, "Gene"],
        how="left",
    )

    # Rename and clean up
    enhanced_df = enhanced_df.rename(columns={"expression_norm": "target_expression"})
    enhanced_df = enhanced_df.drop(columns=["Gene"])

    # 3. Get 'rnase_expression' (Specific 'RNASEH1' Merge)
    # ---------------------------------------------------------
    # Filter master list for RNASEH1
    rnase_vals = expression_master[expression_master["Gene"] == "RNASEH1"]

    # Merge only on cell line ID
    enhanced_df = enhanced_df.merge(
        rnase_vals[[CELL_LINE_DEPMAP, "expression_norm"]],
        on=CELL_LINE_DEPMAP,
        how="left",
        suffixes=("", "_rnase"),  # Safety against collisions
    )

    # Rename
    enhanced_df = enhanced_df.rename(columns={"expression_norm": "rnase_expression"})

    # Define the list of features added
    new_features = ["target_expression", "rnase_expression"]

    return enhanced_df, new_features


def populate_transfection(data):
    """
    One-hot encodes the specified column into 0s and 1s, ensures all expected
    features exist (even if absent from the data), and merges them.
    """
    # The master list of expected categories
    features = ["Electroporation", "Gymnosis", "Lipofection", "Other"]

    # 1. Create binary columns as 0/1 for whatever actually exists in the data
    binary_features = pd.get_dummies(data["transfection_method"], dtype=int)

    # 2. Force the dataframe to have exactly our expected columns, filling missing ones with 0
    binary_features = binary_features.reindex(columns=features, fill_value=0)

    # 3. Join them back to your original dataframe
    data = pd.concat([data, binary_features], axis=1)

    return data, features
