"""TAUSO-side glue for PFRED.

The tauso-FREE PFRED runner lives in the PFREDIntegration submodule
(``notebooks/competitors/PFREDIntegration/integration/pfred_runner.py``). This module
wires it to TAUSO: load the oligo dataset, score it via the runner, write the
``PFRED_SVM`` / ``PFRED_PLS`` features into the TAUSO feature store, and evaluate per
cohort. The runner functions are re-exported here so callers import everything from one
place (``notebooks.competitors.pfred_glue``).

If the import below fails, the submodule isn't initialized:
    git submodule update --init notebooks/competitors/PFREDIntegration
"""

import sys
from pathlib import Path

import pandas as pd

# --- import the tauso-free runner from the PFREDIntegration submodule ---
_RUNNER_DIR = Path(__file__).resolve().parent / "PFREDIntegration" / "integration"
if not (_RUNNER_DIR / "pfred_runner.py").exists():
    raise ImportError(
        "PFREDIntegration submodule not initialized — the PFRED runner is missing.\n"
        "Run: git submodule update --init notebooks/competitors/PFREDIntegration"
    )
sys.path.insert(0, str(_RUNNER_DIR))
from pfred_runner import populate_pfred, score_sequences, validate_docker_container  # noqa: E402,F401

from notebooks.consts import OLIGO_CSV_PROCESSED  # noqa: E402
from notebooks.features.feature_extraction import save_feature  # noqa: E402
from tauso.data.consts import CANONICAL_GENE, CELL_LINE, INHIBITION, SEQUENCE  # noqa: E402


def populate_and_save_pfred(version="oligo"):
    """Load the oligo dataset, score with PFRED, and save features to the TAUSO store."""
    oligo_data = pd.read_csv(OLIGO_CSV_PROCESSED)
    df, features = populate_pfred(oligo_data, seq_col=SEQUENCE)
    for feature in features:
        save_feature(df, feature, version=version)
    return df, features


def median_spearman_per_cohort(df, feature_cols=("PFRED_SVM", "PFRED_PLS"), target_col=INHIBITION):
    """Median per-cohort (gene x cell line) Spearman correlation of each feature vs target."""
    cohort = df[CANONICAL_GENE] + "_" + df[CELL_LINE]
    results = {}
    for col in feature_cols:
        corr = df.groupby(cohort).apply(lambda x: x[col].corr(x[target_col], method="spearman"))
        results[col] = corr.median()
    return pd.Series(results)
