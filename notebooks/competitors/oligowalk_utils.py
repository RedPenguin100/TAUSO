import os
import subprocess
import tempfile
import multiprocessing
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from Bio.Seq import Seq
from tqdm import tqdm

# --- User Imports ---
from tauso.data.consts import SEQUENCE, CANONICAL_GENE

# =====================================================================
# MULTIPROCESSING GLOBALS
# =====================================================================
_WORKER_GENE_DATA = None
_WINDOW_SIZE = 800
_DATAPATH = ""


def _init_worker(gene_data, window_size, datapath):
    """Initializes global variables once per worker process."""
    global _WORKER_GENE_DATA, _WINDOW_SIZE, _DATAPATH
    _WORKER_GENE_DATA = gene_data
    _WINDOW_SIZE = window_size
    _DATAPATH = datapath


# =====================================================================
# HELPER FUNCTIONS
# =====================================================================


def _infer_datapath():
    """Infers RNAstructure DATAPATH from environment or Conda prefix."""
    if "DATAPATH" in os.environ:
        return os.environ["DATAPATH"]

    # Check active conda environment
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        candidate = os.path.join(conda_prefix, "share", "rnastructure", "data_tables")
        # RNAstructure often expects a trailing slash
        if os.path.isdir(candidate):
            return candidate + "/"

    raise EnvironmentError("Could not infer DATAPATH. Please set the DATAPATH environment variable.")


def get_target_site_seq(aso_seq_str):
    """Converts ASO to Target Site (RevComp)."""
    clean_seq = "".join(c for c in str(aso_seq_str) if c.isalpha())
    return str(Seq(clean_seq).reverse_complement())


def _run_single_aso_worker(row_tuple):
    """Core logic to run OligoWalk for one row."""
    idx, row = row_tuple

    metrics = {
        "OW_Overall": None,
        "OW_Duplex": None,
        "OW_Tm": None,
        "OW_Break_Target": None,
        "OW_Intra_Oligo": None,
        "error": None,
    }

    aso_seq = row.get(SEQUENCE)
    gene_name = row.get(CANONICAL_GENE)

    # 1. Validation
    if not gene_name or gene_name not in _WORKER_GENE_DATA:
        metrics["error"] = "Gene missing"
        return idx, metrics

    try:
        full_target_seq = _WORKER_GENE_DATA[gene_name].full_mrna
    except AttributeError:
        metrics["error"] = "No mRNA seq"
        return idx, metrics

    target_site = get_target_site_seq(aso_seq)

    # 2. Find Window
    target_idx = full_target_seq.upper().find(target_site.upper())
    if target_idx == -1:
        metrics["error"] = "Site not found"
        return idx, metrics

    start = max(0, target_idx - _WINDOW_SIZE // 2)
    end = min(len(full_target_seq), target_idx + len(target_site) + _WINDOW_SIZE // 2)
    window_seq = full_target_seq[start:end]
    expected_pos = (target_idx - start) + 1

    # 3. Prepare Files (Use PID to avoid worker collisions)
    pid = os.getpid()
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=f"_{pid}.fasta") as tmp_in:
        tmp_in.write(f">window\n{window_seq}\n")
        tmp_in_name = tmp_in.name
    tmp_out_name = tmp_in_name + ".out"

    # 4. Run Tool
    cmd = ["OligoWalk", tmp_in_name, tmp_out_name, "-m", "1", "-l", str(len(aso_seq))]

    my_env = os.environ.copy()
    my_env["DATAPATH"] = _DATAPATH

    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, env=my_env)

        # 5. Parse Output
        if os.path.exists(tmp_out_name):
            with open(tmp_out_name, "r") as f:
                for line in f:
                    parts = line.split()
                    if len(parts) < 7 or not parts[0].isdigit():
                        continue

                    if int(parts[0]) == expected_pos:
                        metrics["OW_Overall"] = float(parts[2])
                        metrics["OW_Duplex"] = float(parts[3])
                        metrics["OW_Tm"] = float(parts[4])
                        metrics["OW_Break_Target"] = float(parts[5])
                        metrics["OW_Intra_Oligo"] = float(parts[6])
                        break
        else:
            metrics["error"] = "Output file missing"

        if metrics["OW_Overall"] is None and metrics["error"] is None:
            metrics["error"] = f"Pos {expected_pos} not in report"

    except Exception as e:
        metrics["error"] = str(e)
    finally:
        # 6. Cleanup
        if os.path.exists(tmp_in_name):
            os.remove(tmp_in_name)
        if os.path.exists(tmp_out_name):
            os.remove(tmp_out_name)

    return idx, metrics


def populate_oligowalk(df: pd.DataFrame, gene_data: dict, window_size: int = 800, num_workers: int = None):
    """
    Runs OligoWalk sequentially across a DataFrame and appends the features.

    Args:
        df (pd.DataFrame): Input dataframe containing sequences and genes.
        gene_data (dict): Mapping of gene names to gene objects.
        window_size (int): Size of the mRNA window to query.
        num_workers (int): Cores to use (Defaults to cpu_count - 2).

    Returns:
        tuple: (updated_dataframe, list_of_new_features)
    """
    new_features = ["OW_Overall", "OW_Duplex", "OW_Tm", "OW_Break_Target", "OW_Intra_Oligo", "error"]

    if df.empty:
        return df, new_features

    if num_workers is None:
        num_workers = max(1, multiprocessing.cpu_count() - 2)

    datapath = _infer_datapath()
    print(f"Inferred DATAPATH: {datapath}")
    print(f"Running OligoWalk on {len(df)} candidates using {num_workers} parallel workers...")

    # Convert dataframe rows to list of tuples (idx, dict)
    rows_to_process = [(idx, row.to_dict()) for idx, row in df.iterrows()]
    results_dict = {}

    # Execute Parallel Pool
    with ProcessPoolExecutor(
        max_workers=num_workers, initializer=_init_worker, initargs=(gene_data, window_size, datapath)
    ) as executor:
        for idx, metrics in tqdm(
            executor.map(_run_single_aso_worker, rows_to_process), total=len(rows_to_process), desc="OligoWalk Progress"
        ):
            results_dict[idx] = metrics

    # Build results dataframe and join back
    results_df = pd.DataFrame.from_dict(results_dict, orient='index')
    updated_df = df.join(results_df)

    return updated_df, new_features