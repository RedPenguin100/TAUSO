import logging
import os

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Files sourced from ATtRACT.zip (https://attract.cnic.es/attract/static/ATtRACT.zip)
_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


def load_attract_data():
    """
    Parses the bundled ATtRACT files.
    """

    # Point to the bundled files
    csv_path = os.path.join(_DATA_DIR, "RBS_motifs_Homo_sapiens.csv")
    pwm_path = os.path.join(_DATA_DIR, "pwm.txt")

    # Check if they exist before running
    if not os.path.exists(csv_path) or not os.path.exists(pwm_path):
        raise FileNotFoundError(f"Bundled ATtRACT files not found in {_DATA_DIR}.")

    logger.info("Loading RBP metadata from %s...", csv_path)

    # NOTE: We now read the CSV created by setup_attract (comma-separated)
    df = pd.read_csv(csv_path)

    # Map Gene Name -> List of Matrix IDs
    # (No need to filter for Organism here, it was done in setup_attract)
    rbp_to_matrices = df.groupby("Gene_name", observed=True)["Matrix_id"].apply(list).to_dict()

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

    logger.info("Loaded %d PWMs covering %d unique RBPs.", len(pwms), len(rbp_to_matrices))
    return rbp_to_matrices, pwms
