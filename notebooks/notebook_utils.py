
import numpy as np

from tauso.data.consts import CANONICAL_GENE_NAME, INHIBITION_PERCENT


def log_correction(df, correction=0.01):
    df.loc[:, 'log_inhibition'] = -np.log(-df[INHIBITION_PERCENT] + 100 + correction)  # to avoid log 0


def get_unique_genes(df):
    """
    Extracts unique genes from a dataframe, excluding blacklist items.
    """
    genes = df[CANONICAL_GENE_NAME].unique()
    genes_to_remove = {"HBV", "negative_control", np.nan}

    # Filter using list comprehension
    genes_clean = [g for g in genes if g not in genes_to_remove]
    return genes_clean
