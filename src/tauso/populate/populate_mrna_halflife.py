import logging

import pandas as pd

from ..data.consts import CANONICAL_GENE, CELL_LINE
from ..features.context.mrna_halflife import cell_line_mapping

logger = logging.getLogger(__name__)


def populate_mrna_halflife_features(all_data, provider):
    """
    Enriches the dataframe with mRNA half-life features by mapping cell lines
    to TTDB proxies and querying the HalfLifeProvider.

    Args:
        all_data (pd.DataFrame): Dataframe containing CANONICAL_GENE and CELL_LINE columns.
                      ['mRNA_HalfLife', 'HalfLife_Source', 'Mapped_Cell_Proxy']
    """
    logger.info("Calculating stability features for %d rows...", len(all_data))

    features = ["mRNA_HalfLife", "HalfLife_Source", "Mapped_Cell_Proxy"]

    # 2. Define the logic for a single row
    def _get_halflife_features(row):
        # Extract keys using global constants
        gene = row[CANONICAL_GENE]
        cell = row[CELL_LINE]

        # A. Map to TTDB Proxy
        # Uses the imported 'cell_line_mapping' dict to standardize names
        proxy_cell = cell_line_mapping.get(cell, cell)

        # B. Query the Provider
        # Returns: (half_life, source, n_support, std_dev)
        hl_val, source, n, std = provider.get_halflife(gene, proxy_cell)

        # C. Clip Artifacts (Standardize max duration to 48h)
        hl_final = min(hl_val, 48.0)

        # Return the columns to be added
        return pd.Series([hl_final, source, proxy_cell], index=features)

    # 3. Apply to Main DataFrame
    # axis=1 passes the row to the function
    features_df = all_data.apply(_get_halflife_features, axis=1)

    # 4. Merge and Return
    # Concatenate the new columns to the main dataframe matching indices
    enriched_df = pd.concat([all_data, features_df], axis=1)

    logger.info("Features populated successfully.")
    return enriched_df, features
