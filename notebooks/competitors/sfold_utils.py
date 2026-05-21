import os
import shutil
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any, Dict, Optional, Tuple

import numpy as np
import pandas as pd
from Bio.Seq import Seq
from tqdm import tqdm

from tauso.data.consts import CANONICAL_GENE, SEQUENCE

# --- CONFIGURATION & CONSTANTS ---
CONTEXT_WINDOW = 50
SFOLD_PROB_COLUMN_INDEX = 3
SFOLD_BINARY = "sfold"  # Ensure 'sfold' is in your PATH


def process_single_aso(aso_seq: str, full_mrna: str) -> Tuple[float, str]:
    """
    Returns:
        (avg_score, status_message):
        - If success: (score, "Success")
        - If failure: (NaN, "Error: <reason>")
    """
    try:
        # 1. Normalize Sequences (Crucial for T/U and case mismatch)
        norm_mrna = full_mrna.upper().replace("U", "T")
        norm_aso = aso_seq.upper().replace("U", "T")

        target_binding_site = str(Seq(norm_aso).reverse_complement())
        start_pos_global = norm_mrna.find(target_binding_site)

        # --- DIAGNOSTIC DEBUG BLOCK ---
        if start_pos_global == -1:
            # Check if it exists as an exact forward match (sense strand)
            if norm_mrna.find(norm_aso) != -1:
                return np.nan, "Error: Found EXACT match (ASO appears to be sense, not antisense)"

            # Check if it exists as a strict reverse (not complemented)
            if norm_mrna.find(norm_aso[::-1]) != -1:
                return np.nan, "Error: Found REVERSED match (not complemented)"

            # If truly not found, return a snippet to verify mRNA data integrity
            mrna_snippet = norm_mrna[:20] + "..." if len(norm_mrna) > 20 else norm_mrna
            debug_msg = f"Error: Site {target_binding_site} missing. mRNA starts with: {mrna_snippet}"
            return np.nan, debug_msg
        # ------------------------------

        end_pos_global = start_pos_global + len(target_binding_site)

        # 2. Extract Context Window
        # Using original full_mrna for Sfold to preserve intended alphabet if it expects RNA
        w_start = max(0, start_pos_global - CONTEXT_WINDOW)
        w_end = min(len(full_mrna), end_pos_global + CONTEXT_WINDOW)
        subsequence = full_mrna[w_start:w_end]

        # 3. Run Sfold
        temp_dir_base = "/dev/shm" if os.path.exists("/dev/shm") else None

        with tempfile.TemporaryDirectory(dir=temp_dir_base) as tmpdirname:
            input_path = os.path.join(tmpdirname, "input.fa")
            with open(input_path, "w") as f:
                f.write(f">win\n{subsequence}\n")

            cmd = f"{SFOLD_BINARY} -o {tmpdirname} {input_path}"

            try:
                subprocess.run(
                    cmd, shell=True, check=True, cwd=tmpdirname, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
                )
            except subprocess.CalledProcessError as e:
                return np.nan, f"Error: sfold execution failed (Return Code: {e.returncode})"

            # 4. Parse Output
            target_file = os.path.join(tmpdirname, "sstrand.out")
            if not os.path.exists(target_file):
                return np.nan, "Error: sstrand.out file not generated"

            probs = []
            with open(target_file, "r") as f:
                for line in f:
                    parts = line.split()
                    if len(parts) > 4 and parts[0].isdigit():
                        try:
                            probs.append(float(parts[SFOLD_PROB_COLUMN_INDEX]))
                        except ValueError:
                            continue

            if not probs:
                return np.nan, "Error: Parsed probability list is empty"

            profile = np.array(probs)

        # 5. Calculate Score
        local_start = start_pos_global - w_start
        local_end = local_start + len(target_binding_site)

        if local_end <= len(profile):
            region_probs = profile[local_start:local_end]
            return np.mean(region_probs), "Success"
        else:
            return np.nan, f"Error: Index out of bounds (Profile len: {len(profile)}, req end: {local_end})"

    except Exception as e:
        return np.nan, f"Error: Uncaught exception ({str(e)})"


# --- MAIN WRAPPER FUNCTION ---
def calculate_sfold_accessibility(
    df: pd.DataFrame,
    gene_to_data: Dict[str, Any],
    score_col: str = "sfold_accessibility",
    error_col: str = "sfold_error",
    n_cores: Optional[int] = None,
) -> pd.DataFrame:
    if n_cores is None:
        n_cores = max(1, os.cpu_count() - 2)

    if shutil.which(SFOLD_BINARY) is None:
        raise RuntimeError(f"The '{SFOLD_BINARY}' executable was not found in PATH.")

    tasks = []
    unique_genes = df[CANONICAL_GENE].unique()

    for gene in unique_genes:
        if gene not in gene_to_data:
            continue

        gene_obj = gene_to_data[gene]
        full_seq = gene_obj.full_mrna if hasattr(gene_obj, "full_mrna") else str(gene_obj)

        if not full_seq or len(full_seq) < 10:
            continue

        gene_mask = df[CANONICAL_GENE] == gene
        for idx, row in df[gene_mask].iterrows():
            aso_seq = row[SEQUENCE].strip()
            tasks.append((idx, aso_seq, full_seq))

    print(f"Queueing {len(tasks)} jobs on {n_cores} cores...")

    if score_col not in df.columns:
        df[score_col] = np.nan
    df[error_col] = "Pending"

    results_map = {}

    with ProcessPoolExecutor(max_workers=n_cores) as executor:
        future_to_idx = {executor.submit(process_single_aso, t[1], t[2]): t[0] for t in tasks}

        for future in tqdm(as_completed(future_to_idx), total=len(tasks), desc="Sfold Analysis"):
            idx = future_to_idx[future]
            try:
                score, message = future.result()
                results_map[idx] = (score, message)
            except Exception as e:
                results_map[idx] = (np.nan, f"Error: Critical wrapper failure ({str(e)})")

    print("Updating DataFrame...")

    for idx, (score, msg) in results_map.items():
        df.at[idx, score_col] = score
        df.at[idx, error_col] = msg

    return df
