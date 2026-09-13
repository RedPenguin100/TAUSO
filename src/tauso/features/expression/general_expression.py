import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

logger = logging.getLogger(__name__)

GENE_COLUMN_BATCH = 1000
"""Gene columns averaged at once. The table is ~60,000 columns wide and only the per-column
mean is kept, so holding a slice at a time costs a fraction of holding the table."""


def get_general_expression_of_genes(EXP_path, valid_genes):
    """
    Loads expression data, filters for valid genes based on GTF (e.g., protein_coding),
    averages across cell lines, and converts to TPM.
    Expects a Parquet file (converted from the original DepMap CSV by setup-depmap).
    """
    EXP_path = Path(EXP_path)
    parquet_path = EXP_path.with_suffix(".parquet")

    valid_genes_set = set(valid_genes)

    handle = pq.ParquetFile(parquet_path)
    all_cols = handle.schema.names
    model_col = "ModelID" if "ModelID" in all_cols else all_cols[0]

    potential_gene_cols = [c for c in all_cols if "(" in c and c != model_col]
    valid_cols = [c for c in potential_gene_cols if c.split(" (")[0] in valid_genes_set]

    logger.info(f"Filtering: Kept {len(valid_cols)} genes out of {len(potential_gene_cols)} total columns.")

    if not valid_cols:
        return pd.DataFrame(columns=["Gene", "expression_norm", "expression_TPM"])

    mean_exp = np.empty(len(valid_cols))
    for start in range(0, len(valid_cols), GENE_COLUMN_BATCH):
        batch = handle.read(columns=valid_cols[start : start + GENE_COLUMN_BATCH])
        for offset, column in enumerate(batch.columns):
            mean_exp[start + offset] = np.nanmean(column.to_numpy(zero_copy_only=False))

    mean_exp_data = pd.DataFrame({"Gene": valid_cols, "expression_norm": mean_exp})
    mean_exp_data["expression_TPM"] = 2 ** mean_exp_data["expression_norm"]

    return mean_exp_data.sort_values("expression_norm", ascending=False)
