"""Shared helpers for the competitor-score regeneration scripts."""
import pandas as pd
from scipy.stats import spearmanr

from notebooks.consts import OLIGO_CSV_INDEXED
from notebooks.notebook_utils import get_unique_genes
from notebooks.preprocessing import process_oligo_data_rename
from tauso.data.consts import CANONICAL_GENE_NAME, CELL_LINE, INHIBITION_PERCENT
from tauso.genome.read_human_genome import get_locus_to_data_dict


def load_indexed():
    """Indexed table standardized to canonical names, restricted to gene-mapped rows."""
    df = process_oligo_data_rename(pd.read_csv(OLIGO_CSV_INDEXED))
    return df[df[CANONICAL_GENE_NAME].notna()].reset_index(drop=True)


def load_gene_mapped():
    """load_indexed() plus the gene->mRNA dict (for tools that need the target sequence)."""
    df = load_indexed()
    gene_to_data = get_locus_to_data_dict(include_introns=True, gene_subset=get_unique_genes(df))
    return df, gene_to_data


def validate_and_report(df, feature_cols):
    """Fail loudly if a feature came out degenerate (tool failure); else report per-cohort
    (gene x cell) median Spearman vs inhibition so a bad run is obvious."""
    for col in feature_cols:
        if df[col].nunique(dropna=True) <= 1:
            raise SystemExit(f"{col}: degenerate output (all NaN/constant) -- the tool likely failed "
                             f"(check binary / version / DATAPATH).")

    if INHIBITION_PERCENT not in df.columns:
        return
    cohort = df[CANONICAL_GENE_NAME].astype(str) + "_" + df[CELL_LINE].astype(str)
    for col in feature_cols:
        def cohort_spearman(g, _col=col):
            g = g[[_col, INHIBITION_PERCENT]].dropna()
            return spearmanr(g[_col], g[INHIBITION_PERCENT]).correlation if len(g) >= 5 else float("nan")
        med = df.groupby(cohort, group_keys=False).apply(cohort_spearman).median()
        print(f"  {col}: per-cohort median Spearman vs inhibition = {med:.4f}")
