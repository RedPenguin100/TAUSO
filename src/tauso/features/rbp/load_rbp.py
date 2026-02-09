import numpy as np
import pandas as pd
import os

from ...data.data import get_data_dir


def load_attract_data():
    """
    Parses the processed ATtRACT files.
    """

    data_dir = get_data_dir()

    # Point to the files created by the CLI
    csv_path = os.path.join(data_dir, "RBS_motifs_Homo_sapiens.csv")
    pwm_path = os.path.join(data_dir, "pwm.txt")

    # Check if they exist before running
    if not os.path.exists(csv_path) or not os.path.exists(pwm_path):
        raise FileNotFoundError("ATtRACT files not found. Please run 'tauso setup-attract' first.")

    print(f"Loading RBP metadata from {csv_path}...")

    # NOTE: We now read the CSV created by setup_attract (comma-separated)
    df = pd.read_csv(csv_path)

    # Map Gene Name -> List of Matrix IDs
    # (No need to filter for Organism here, it was done in setup_attract)
    rbp_to_matrices = df.groupby('Gene_name')['Matrix_id'].apply(list).to_dict()

    print(f"Loading PWM matrices from {pwm_path}...")
    pwms = {}
    current_matrix_id = None
    matrix_data = []

    with open(pwm_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue

            if line.startswith('>'):
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

    print(f"Loaded {len(pwms)} PWMs covering {len(rbp_to_matrices)} unique RBPs.")
    return rbp_to_matrices, pwms
