"""Shared helpers for the competitor-score regeneration scripts."""
import pandas as pd

from notebooks.consts import OLIGO_CSV_INDEXED
from notebooks.notebook_utils import get_unique_genes
from notebooks.preprocessing import process_oligo_data_rename
from tauso.data.consts import CANONICAL_GENE_NAME
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


def validate_features(df, feature_cols):
    """Fail loudly if a feature came out degenerate (all NaN/constant) -- the tool likely failed."""
    for col in feature_cols:
        if df[col].nunique(dropna=True) <= 1:
            raise SystemExit(f"{col}: degenerate output (all NaN/constant) -- the tool likely failed "
                             f"(check binary / version / DATAPATH).")
