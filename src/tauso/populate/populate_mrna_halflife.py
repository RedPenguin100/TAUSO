import logging

import numpy as np
import pandas as pd

from ..data.consts import CANONICAL_GENE_NAME, CELL_LINE
from ..features.context.ttdb_cell_names import cell_name_to_ttdb

logger = logging.getLogger(__name__)


def _halflife_for(provider, gene, cell):
    """One gene in one cell line: the clipped half-life, its provenance, and the cell used."""
    # Fix the cell name to TTDB's spelling; unknown names pass through and resolve to gene-level.
    ttdb_cell = cell_name_to_ttdb(cell)

    res = provider.get_halflife(gene, ttdb_cell)

    # Standardize max duration to 48h, preserving NaN for absent genes.
    hl_final = min(res.half_life, 48.0) if np.isfinite(res.half_life) else np.nan

    return hl_final, res.source, ttdb_cell


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

    # The answer depends on nothing but the gene and the cell line, and the corpus holds far
    # fewer distinct pairs than rows, so each pair is looked up once and mapped back.
    pairs = all_data[[CANONICAL_GENE_NAME, CELL_LINE]].drop_duplicates()
    pairs[features] = pd.DataFrame(
        [_halflife_for(provider, gene, cell) for gene, cell in pairs.itertuples(index=False)],
        columns=features,
        index=pairs.index,
    )

    per_row = all_data[[CANONICAL_GENE_NAME, CELL_LINE]].merge(pairs, on=[CANONICAL_GENE_NAME, CELL_LINE], how="left")
    for feature in features:
        all_data[feature] = per_row[feature].to_numpy()

    logger.info("Features populated successfully.")
    return all_data, features
