import logging
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

from ..data.consts import CANONICAL_GENE_NAME
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

