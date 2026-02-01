import numpy as np
import pandas as pd

from ..data.consts import CANONICAL_GENE
from ..features.codon_usage.tai import calc_tAI, tai_weights

def populate_tai(df: pd.DataFrame, cds_windows: list, registry: dict) -> tuple[pd.DataFrame, list[str]]:
    """
    Calculates local and global tAI scores.
    Raises KeyError if required local context columns are missing.
    """
    weights = tai_weights("hm")
    feature_names = []

    # 1. Local tAI Scores
    print("Calculating local tAI scores for CDS windows...")
    for flank in cds_windows:
        local_col = f"local_coding_region_around_ASO_{flank}"

        # Strict check: Raise error if column is missing
        if local_col not in df.columns:
            raise KeyError(f"Missing required column for tAI calculation: '{local_col}'. Please run add_external_mrna_and_context_columns")

        feat_name = f"tAI_score_{flank}_CDS"
        df[feat_name] = df[local_col].apply(lambda seq: calc_tAI(seq, weights))
        feature_names.append(feat_name)

    # 2. Global tAI Scores (Registry Lookup)
    print(f"Pre-calculating Global tAI for {len(registry)} unique genes...")
    gene_tai_lookup = {
        gene: calc_tAI(data.get('cds_sequence'), weights) if data.get('cds_sequence') else np.nan
        for gene, data in registry.items()
    }

    print("Mapping Global tAI scores to dataframe...")
    df["tAI_score_global_CDS"] = df[CANONICAL_GENE].map(gene_tai_lookup)
    feature_names.append("tAI_score_global_CDS")

    return df, feature_names