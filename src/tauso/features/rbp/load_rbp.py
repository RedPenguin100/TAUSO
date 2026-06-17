import logging
import os

import numpy as np
import pandas as pd

from ...data.data import get_data_dir

logger = logging.getLogger(__name__)

# Files sourced from a frozen ATtRACT mirror on Zenodo (https://zenodo.org/records/20366079).
# Downloaded into TAUSO_DATA_DIR via `tauso setup-attract`.
_DATA_DIR = os.path.join(get_data_dir(), "attract")

_MIN_MOTIF_LEN = 6
# Minimum total information content (bits) a PWM must have to be kept. Increase to restrict the
# RBP set to higher-confidence motifs; 0 keeps all of them.
_MIN_TOTAL_IC = 0.0
_DB_RANK = {"S": 3, "R": 2, "C": 1}


def _pwm_total_ic(pwm):
    """Total information content (bits) of a PWM against a uniform background."""
    eps = 1e-9
    p = pwm + eps
    p = p / p.sum(axis=1, keepdims=True)
    return float(np.sum(p * np.log2(p / 0.25), axis=1).sum())


def _rank_matrices_per_gene(df, pwms):
    """Map Gene_name -> matrix IDs ordered best-first.

    Filters out PWMs below `_MIN_TOTAL_IC` bits, then ranks survivors by
    (total_IC, length >= _MIN_MOTIF_LEN, Score, DB tier, length).
    Genes whose best matrix is below the IC threshold are dropped entirely.
    """
    score = pd.to_numeric(df["Score"].astype(str).str.replace("*", "", regex=False), errors="coerce")
    q = (
        pd.DataFrame({"gene": df["Gene_name"], "matrix_id": df["Matrix_id"], "score": score, "db": df["Database"]})
        .groupby(["gene", "matrix_id"], observed=True)
        .agg(score=("score", "max"), db=("db", "first"))
        .reset_index()
    )
    pwm_ic = {matrix_id: _pwm_total_ic(pwm) for matrix_id, pwm in pwms.items()}
    rbp_to_matrices = {}
    for gene, sub in q.groupby("gene", observed=True):
        ranked = []
        for _, row in sub.iterrows():
            matrix_id = row["matrix_id"]
            if matrix_id not in pwms:
                continue
            ic = pwm_ic[matrix_id]
            if ic < _MIN_TOTAL_IC:
                continue
            length = pwms[matrix_id].shape[0]
            sc = -1.0 if pd.isna(row["score"]) else float(row["score"])
            ranked.append((ic, length >= _MIN_MOTIF_LEN, sc, _DB_RANK.get(row["db"], 0), length, matrix_id))
        if not ranked:
            continue
        ranked.sort(reverse=True)
        rbp_to_matrices[gene] = [r[-1] for r in ranked]
    return rbp_to_matrices


def load_attract_data():
    """
    Parses the bundled ATtRACT files.
    """

    # Point to the bundled files
    csv_path = os.path.join(_DATA_DIR, "RBS_motifs_Homo_sapiens.csv")
    pwm_path = os.path.join(_DATA_DIR, "pwm.txt")

    # Check if they exist before running
    if not os.path.exists(csv_path) or not os.path.exists(pwm_path):
        raise FileNotFoundError(f"ATtRACT files not found in {_DATA_DIR}. Run 'tauso setup-attract' to download them.")

    logger.info("Loading RBP metadata from %s...", csv_path)
    df = pd.read_csv(csv_path)

    logger.info("Loading PWM matrices from %s...", pwm_path)
    pwms = {}
    current_matrix_id = None
    matrix_data = []

    with open(pwm_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_matrix_id and matrix_data:
                    pwms[current_matrix_id] = np.array(matrix_data)

                # Header format: >MatrixID Length
                parts = line.split()
                current_matrix_id = parts[0][1:]  # Remove '>'
                matrix_data = []
            else:
                # Cols: A C G U
                probs = [float(x) for x in line.split()]
                matrix_data.append(probs)

        # Save the last matrix
        if current_matrix_id and matrix_data:
            pwms[current_matrix_id] = np.array(matrix_data)

    rbp_to_matrices = _rank_matrices_per_gene(df, pwms)

    logger.info("Loaded %d PWMs covering %d unique RBPs.", len(pwms), len(rbp_to_matrices))
    return rbp_to_matrices, pwms
