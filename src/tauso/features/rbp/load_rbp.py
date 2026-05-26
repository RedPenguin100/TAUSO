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
_DB_RANK = {"S": 3, "R": 2, "C": 1}


def _rank_matrices_per_gene(df, pwms):
    """Map Gene_name -> matrix IDs ordered best-first (Score, length >= _MIN_MOTIF_LEN, assay)."""
    score = pd.to_numeric(df["Score"].astype(str).str.replace("*", "", regex=False), errors="coerce")
    q = (
        pd.DataFrame({"gene": df["Gene_name"], "mid": df["Matrix_id"], "score": score, "db": df["Database"]})
        .groupby(["gene", "mid"], observed=True)
        .agg(score=("score", "max"), db=("db", "first"))
        .reset_index()
    )
    rbp_to_matrices = {}
    for gene, sub in q.groupby("gene", observed=True):
        ranked = []
        for _, row in sub.iterrows():
            mid = row["mid"]
            if mid not in pwms:
                continue
            length = pwms[mid].shape[0]
            sc = -1.0 if pd.isna(row["score"]) else float(row["score"])
            ranked.append((length >= _MIN_MOTIF_LEN, sc, _DB_RANK.get(row["db"], 0), length, mid))
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

    # NOTE: We now read the CSV created by setup_attract (comma-separated)
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
