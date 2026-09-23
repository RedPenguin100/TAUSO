import logging
from pathlib import Path

import pandas as pd

from ..frames import ARROW_STRINGS, release_arrow_memory
from .depmap_parquet import expression_columns, mean_expression, parquet_exists

logger = logging.getLogger(__name__)

# The DepMap matrix every gene's mean is taken over, and the table of means itself. The
# means depend on that file and on which genes the annotation calls valid, so they are
# computed once by `tauso build-general-expression` rather than per run.
OMICS_FILENAME = "OmicsExpressionTPMLogp1HumanAllGenesStranded.parquet"
GENERAL_EXPRESSION_FILENAME = "general_expression.parquet"


def get_general_expression_of_genes(EXP_path, valid_genes):
    """
    Loads expression data, filters for valid genes based on GTF (e.g., protein_coding),
    averages across cell lines, and converts to TPM.
    Expects the Parquet setup-depmap converted the original DepMap CSV to; EXP_path names the
    CSV or the Parquet, and the table is found from it.
    """
    valid_genes_set = set(valid_genes)

    if not parquet_exists(EXP_path):
        raise FileNotFoundError(f"No Parquet for {EXP_path}. Run 'tauso setup-depmap'.")
    potential_gene_cols = [c for c in expression_columns(EXP_path) if "(" in c]
    valid_cols = [c for c in potential_gene_cols if c.split(" (")[0] in valid_genes_set]

    logger.info(f"Filtering: Kept {len(valid_cols)} genes out of {len(potential_gene_cols)} total columns.")

    if not valid_cols:
        return pd.DataFrame(columns=["Gene", "expression_norm", "expression_TPM"])

    mean_exp = mean_expression(EXP_path, valid_cols)

    mean_exp_data = pd.DataFrame({"Gene": pd.array(valid_cols, dtype=ARROW_STRINGS), "expression_norm": mean_exp})
    mean_exp_data["expression_TPM"] = 2 ** mean_exp_data["expression_norm"]

    return mean_exp_data.sort_values("expression_norm", ascending=False)


def general_expression_path() -> Path:
    from ..data.data import get_data_dir

    return Path(get_data_dir()) / GENERAL_EXPRESSION_FILENAME


def build_general_expression(genome: str = "GRCh38") -> pd.DataFrame:
    """Compute the per-gene mean expression across the cohort and write it to the data dir.

    Reading the DepMap matrix costs far more than the table it produces -- 60,655 columns
    against 38,584 rows of means -- and nothing about it changes from run to run.
    """
    import os

    from ..common.gtf import filter_gtf_genes
    from ..data.data import get_data_dir, load_gtf_db

    valid_genes = filter_gtf_genes(load_gtf_db(genome), filter_mode="non_mt")
    table = get_general_expression_of_genes(Path(get_data_dir()) / OMICS_FILENAME, valid_genes)

    # write-then-rename: two runs starting at once must not read a half-written table
    path = general_expression_path()
    tmp = f"{path}.{os.getpid()}.tmp"
    table.to_parquet(tmp, index=False)
    os.replace(tmp, path)
    release_arrow_memory()
    return table


def load_general_expression(genome: str = "GRCh38") -> pd.DataFrame:
    """The per-gene mean expression, most expressed first.

    Built on first use and read from the data dir after; `tauso build-general-expression`
    does the same thing ahead of time.
    """
    from ..frames import arrow_strings

    path = general_expression_path()
    if not path.exists():
        logger.info("No general expression table yet; building %s (once).", path)
        return build_general_expression(genome)
    return arrow_strings(pd.read_parquet(path))
