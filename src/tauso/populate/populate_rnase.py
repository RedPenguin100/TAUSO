import pandas as pd

from ..data.consts import CHEMICAL_PATTERN
from ..features.rnase_motifs.calculate_rnase import (
    add_rnaseh1_motif_features,
    add_rnaseh1_scores_dinuc,
    add_rnaseh1_scores_krel_dinuc,
    add_rnaseh1_scores_krel_nt,
    add_rnaseh1_scores_nt,
)


def populate_rnase_features(df: pd.DataFrame):
    if CHEMICAL_PATTERN not in df.columns:
        raise KeyError(f"Required column '{CHEMICAL_PATTERN}' is missing from the data.")
    df, features_dinuc_krel = add_rnaseh1_scores_krel_dinuc(df)
    df, features_dinuc = add_rnaseh1_scores_dinuc(df)
    df, features_single_krel = add_rnaseh1_scores_krel_nt(df)
    df, features_single = add_rnaseh1_scores_nt(df)
    df, features_motifs = add_rnaseh1_motif_features(df)
    return (df, features_dinuc_krel + features_dinuc + features_single_krel + features_single + features_motifs)
