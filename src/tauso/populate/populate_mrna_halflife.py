import logging

import numpy as np
import pandas as pd

from ..data.consts import CANONICAL_GENE_NAME, CELL_LINE
from ..features.context.ttdb_cell_names import cell_name_to_ttdb

logger = logging.getLogger(__name__)


def populate_mrna_halflife_features(all_data, provider):
    """
    Enriches the dataframe with mRNA half-life features by mapping cell lines
    to TTDB cell types and querying the HalfLifeProvider.

    Adds, per row:
      halflife_value       clipped half-life in hours (NaN if the gene is absent)
      halflife_source      human-readable provenance string
      halflife_cell_proxy  the TTDB cell type the lookup used
    """
    logger.info(f"Calculating stability features for {len(all_data)} rows...")

    features = ["halflife_value", "halflife_source", "halflife_cell_proxy"]

    def _get_halflife_features(row):
        gene = row[CANONICAL_GENE_NAME]
        cell = row[CELL_LINE]

        # Fix the cell name to TTDB's spelling; unknown names pass through and resolve to gene-level.
        ttdb_cell = cell_name_to_ttdb(cell)

        res = provider.get_halflife(gene, ttdb_cell)

        # Standardize max duration to 48h, preserving NaN for absent genes.
        hl_final = min(res.half_life, 48.0) if np.isfinite(res.half_life) else np.nan

        return pd.Series([hl_final, res.source, ttdb_cell], index=features)

    features_df = all_data.apply(_get_halflife_features, axis=1)

    enriched_df = pd.concat([all_data, features_df], axis=1)

    logger.info("Features populated successfully.")
    return enriched_df, features
