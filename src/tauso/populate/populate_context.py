import logging
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from ..data.consts import CANONICAL_GENE_NAME, CELL_LINE_DEPMAP, TRANSFECTION_RAW
from ..features.context.ribo_seq import (
    feature_names,
    get_feature_prefix,
    get_ribo_bigwig_path,
    process_gene_group,
)

logger = logging.getLogger(__name__)


def _run_gene_groups_parallel(bw_path, gene_groups, flanks, how, n_jobs, prefix="ribo"):
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
            return gene, process_gene_group(bw_path, gene_rows, flanks, how, prefix), None
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


def populate_ribo_seq(organism, aso_df, flanks=(0, 10, 20, 50, 100, 125, 150), how="mean", n_jobs=1, track="40s"):
    if organism != "human":
        raise ValueError("Unsupported organism for ribo_seq feature")
    if how not in {"sum", "mean", "max", "nz_mean", "nz_frac"}:
        raise ValueError(how)

    bw_path = str(get_ribo_bigwig_path(track))
    prefix = get_feature_prefix(track)  # "ribo_40s" or "ribo_80s"
    feat_cols = feature_names(flanks, how, prefix=prefix)

    for col in feat_cols:
        aso_df[col] = np.nan

    valid_mask = aso_df["chrom"].notna() & aso_df["target_start"].notna()
    valid_df = aso_df[valid_mask]
    n_valid = len(valid_df)
    n_genes = valid_df[CANONICAL_GENE_NAME].nunique()

    logger.info(
        "populate_ribo_seq[%s]: %d valid rows across %d genes × %d features (n_jobs=%d)",
        track,
        n_valid,
        n_genes,
        len(feat_cols),
        n_jobs,
    )

    gene_groups = list(valid_df.groupby(CANONICAL_GENE_NAME))

    t0 = time.perf_counter()
    results, skipped = _run_gene_groups_parallel(bw_path, gene_groups, flanks, how, n_jobs, prefix=prefix)
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


_RNASE_GENE = "RNASEH1"

# Scavenger receptors that govern gymnotic (naked ASO) uptake — STAB2 and MSR1 (SR-A1)
# are the primary endocytic receptors responsible for gymnosis efficiency across cell lines.
_SCAVENGER_RECEPTOR_GENES = {
    "expr_stab2": "STAB2",
    "expr_msr1": "MSR1",  # SR-A1
}

_SPECIAL_GENES: Dict[str, str] = {"expr_rnase": _RNASE_GENE, **_SCAVENGER_RECEPTOR_GENES}

TARGET_EXPRESSION_FEATURE_NAMES: List[str] = ["expr_target"]
SPECIAL_GENE_EXPRESSION_FEATURE_NAMES: List[str] = list(_SPECIAL_GENES.keys())
EXPRESSION_FEATURE_NAMES: List[str] = TARGET_EXPRESSION_FEATURE_NAMES + SPECIAL_GENE_EXPRESSION_FEATURE_NAMES


def _build_expression_master(expression_dict: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    dfs = []
    for depmap_id, t_df in expression_dict.items():
        temp = t_df[["Gene", "expression_norm"]].copy()
        temp[CELL_LINE_DEPMAP] = depmap_id
        dfs.append(temp)
    master = pd.concat(dfs, ignore_index=True)
    return master.drop_duplicates(subset=[CELL_LINE_DEPMAP, "Gene"])


def _merge_fixed_gene(
    df: pd.DataFrame, expression_master: pd.DataFrame, gene_symbol: str, feat_name: str
) -> pd.DataFrame:
    gene_vals = expression_master[expression_master["Gene"] == gene_symbol]
    merged = df.merge(
        gene_vals[[CELL_LINE_DEPMAP, "expression_norm"]],
        on=CELL_LINE_DEPMAP,
        how="left",
        suffixes=("", f"_{feat_name}"),
    )
    return merged.rename(columns={"expression_norm": feat_name})


def populate_target_expression(
    df: pd.DataFrame, expression_dict: Dict[str, pd.DataFrame]
) -> Tuple[pd.DataFrame, List[str]]:
    """Merges the ASO's target gene expression (per cell line × gene) into df."""
    existing_cols = [c for c in TARGET_EXPRESSION_FEATURE_NAMES if c in df.columns]
    if existing_cols:
        logger.warning("Dropping existing columns to prevent alignment breaks: %s", existing_cols)
        df = df.drop(columns=existing_cols)

    expression_master = _build_expression_master(expression_dict)

    enhanced_df = df.merge(
        expression_master,
        left_on=[CELL_LINE_DEPMAP, CANONICAL_GENE_NAME],
        right_on=[CELL_LINE_DEPMAP, "Gene"],
        how="left",
    )
    enhanced_df = enhanced_df.rename(columns={"expression_norm": "expr_target"})
    enhanced_df = enhanced_df.drop(columns=["Gene"])

    return enhanced_df, TARGET_EXPRESSION_FEATURE_NAMES


def populate_special_gene_expression(
    df: pd.DataFrame, expression_dict: Dict[str, pd.DataFrame]
) -> Tuple[pd.DataFrame, List[str]]:
    """Merges per-cell-line expression for fixed special genes (RNASEH1 + gymnosis scavenger receptors)."""
    existing_cols = [c for c in SPECIAL_GENE_EXPRESSION_FEATURE_NAMES if c in df.columns]
    if existing_cols:
        logger.warning("Dropping existing columns to prevent alignment breaks: %s", existing_cols)
        df = df.drop(columns=existing_cols)

    expression_master = _build_expression_master(expression_dict)

    for feat_name, gene_symbol in _SPECIAL_GENES.items():
        df = _merge_fixed_gene(df, expression_master, gene_symbol, feat_name)

    return df, SPECIAL_GENE_EXPRESSION_FEATURE_NAMES


_TRANSFECTION_LABEL_TO_COLUMN = {
    "Electroporation": "transfection_electroporation",
    "Gymnosis": "transfection_gymnosis",
    "Lipofection": "transfection_lipofection",
}


def populate_transfection(data):
    """
    One-hot encodes the raw transfection_raw label into three columns:
    transfection_electroporation / transfection_gymnosis / transfection_lipofection.
    Rows whose transfection_raw isn't one of those three (e.g. "Other", missing,
    or any unrecognized label) get NaN across all three -- the model treats them
    as missing rather than as confidently "not Electroporation, not Gymnosis,
    not Lipofection".

    Returns float64 columns because NaN cannot live in an int column.
    """
    features = list(_TRANSFECTION_LABEL_TO_COLUMN.values())
    labels = list(_TRANSFECTION_LABEL_TO_COLUMN.keys())

    method = data[TRANSFECTION_RAW]
    binary = pd.DataFrame(
        {col: (method == label).astype("float64") for label, col in _TRANSFECTION_LABEL_TO_COLUMN.items()},
        index=data.index,
    )
    unknown = ~method.isin(labels)
    binary.loc[unknown, features] = np.nan

    data = pd.concat([data, binary], axis=1)
    return data, features
