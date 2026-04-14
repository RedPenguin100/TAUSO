import os
import subprocess
import tempfile
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
from tauso.data.consts import CANONICAL_GENE, SEQUENCE

# --- Configuration ---
MIRANDA_EXEC = "miranda"
MIRANDA_SCORE_THRESHOLD = "50.0"
MIRANDA_ENERGY_THRESHOLD = "0.0"


def _process_gene_group(args):
    """
    Worker function to run miRanda on a single gene group.
    This runs in a separate process.
    """
    gene_id, gene_subset, target_seq = args

    # Use RAM disk if available for zero-latency I/O
    shm_path = "/dev/shm" if os.path.exists("/dev/shm") else None

    results = {}

    with tempfile.TemporaryDirectory(dir=shm_path) as temp_dir:
        # 1. Create Target FASTA
        target_fasta_path = os.path.join(temp_dir, "target.fa")
        with open(target_fasta_path, "w") as f:
            f.write(f">{gene_id}\n{target_seq}\n")

        # 2. Create ASOs FASTA
        aso_fasta_path = os.path.join(temp_dir, "asos.fa")
        with open(aso_fasta_path, "w") as f:
            for idx, row in gene_subset.iterrows():
                if pd.notna(row[SEQUENCE]):
                    f.write(f">{idx}\n{row[SEQUENCE]}\n")

        # 3. Run miRanda
        output_path = os.path.join(temp_dir, "miranda_out.txt")
        cmd = [
            MIRANDA_EXEC,
            aso_fasta_path,
            target_fasta_path,
            "-sc",
            MIRANDA_SCORE_THRESHOLD,
            "-en",
            MIRANDA_ENERGY_THRESHOLD,
            "-out",
            output_path,
            "-strict",
        ]

        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except (subprocess.CalledProcessError, FileNotFoundError):
            return results  # Return empty if it fails

        # 4. Parse Results
        if os.path.exists(output_path):
            with open(output_path, "r") as f:
                for line in f:
                    if line.startswith(">"):
                        parts = line.strip().split()
                        if len(parts) >= 4:
                            try:
                                idx = int(parts[0][1:])
                                score = float(parts[2])
                                energy = float(parts[3])

                                # Keep the Lowest Energy (strongest binding)
                                if idx in results:
                                    curr_score, curr_energy = results[idx]
                                    if energy < curr_energy:
                                        results[idx] = (score, energy)
                                else:
                                    results[idx] = (score, energy)

                            except ValueError:
                                continue

    return results


def get_populated_df_with_miranda_features(df, gene_to_data, max_workers=None):
    """
    Populates the DataFrame with miRanda scores and energy values using multiprocessing.
    """
    df_out = df.copy()
    grouped = df_out.groupby(CANONICAL_GENE)

    # --- Early Error Check ---
    # Fail fast before doing any heavy lifting
    for gene_id, _ in grouped:
        if gene_id not in gene_to_data:
            raise ValueError(f"Gene {gene_id} not in gene_to_data")

    # 1. Prepare tasks for the worker pool
    tasks = []
    print("Preparing data chunks...")
    for gene_id, gene_subset in grouped:
        target_seq = gene_to_data[gene_id].full_mrna

        # Slice only the required column to keep Inter-Process Communication (IPC) fast
        subset_minimal = gene_subset[[SEQUENCE]]
        tasks.append((gene_id, subset_minimal, target_seq))

    print(f"Starting parallel miRanda analysis on {len(df_out)} sequences across {len(tasks)} gene groups...")

    all_results = {}

    # 2. Execute in parallel
    # max_workers=None will automatically use all available CPU cores/threads
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        futures = [executor.submit(_process_gene_group, task) for task in tasks]

        # Track progress as batches complete
        with tqdm(total=len(tasks), desc="Processing Gene Batches") as pbar:
            for future in as_completed(futures):
                try:
                    # Update our master dictionary with the results from the worker
                    worker_results = future.result()
                    all_results.update(worker_results)
                except Exception as e:
                    print(f"Worker failed with error: {e}")

                pbar.update(1)

    # 3. Vectorized DataFrame Update
    print("Mapping results to DataFrame...")

    if all_results:
        # Convert dictionary to a quick DataFrame
        res_df = pd.DataFrame.from_dict(all_results, orient="index", columns=["miranda_score", "miranda_energy"])

        # Map values back based on index, filling missing sequences with 0.0
        df_out["miranda_score"] = df_out.index.map(res_df["miranda_score"]).fillna(0.0)
        df_out["miranda_energy"] = df_out.index.map(res_df["miranda_energy"]).fillna(0.0)
    else:
        # Failsafe if no hits were found across all sequences
        df_out['miranda_score'] = 0.0
        df_out['miranda_energy'] = 0.0

    return df_out