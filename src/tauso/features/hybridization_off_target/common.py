import logging
from pathlib import Path

import pandas as pd

# TODO: move this function to a folder / file about expression

logger = logging.getLogger(__name__)


def get_general_expression_of_genes(EXP_path, valid_genes):
    """
    Loads expression data, filters for valid genes based on GTF (e.g., protein_coding),
    averages across cell lines, and converts to TPM.
    Expects a Parquet file (converted from the original DepMap CSV by setup-depmap).
    """
    import pyarrow.parquet as pq

    EXP_path = Path(EXP_path)
    parquet_path = EXP_path.with_suffix(".parquet")

    valid_genes_set = set(valid_genes)

    # Read schema without loading data to find matching gene columns
    schema = pq.read_schema(parquet_path)
    all_cols = schema.names
    model_col = "ModelID" if "ModelID" in all_cols else all_cols[0]

    potential_gene_cols = [c for c in all_cols if "(" in c and c != model_col]
    valid_cols = [c for c in potential_gene_cols if c.split(" (")[0] in valid_genes_set]

    logger.info("Filtering: Kept %d genes out of %d total columns.", len(valid_cols), len(potential_gene_cols))

    if not valid_cols:
        return pd.DataFrame(columns=["Gene", "expression_norm", "expression_TPM"])

    # Column-selective read — pyarrow only decompresses the requested columns
    exp_data = pd.read_parquet(parquet_path, columns=valid_cols)

    mean_exp = exp_data.mean(axis=0)
    mean_exp_data = mean_exp.reset_index()
    mean_exp_data.columns = ["Gene", "expression_norm"]
    mean_exp_data["expression_TPM"] = 2 ** mean_exp_data["expression_norm"]

    return mean_exp_data.sort_values("expression_norm", ascending=False)
