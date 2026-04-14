import pandas as pd
import subprocess
import os
import re

from tauso.data.consts import SEQUENCE

import subprocess


def validate_docker_container(
    container_name: str = "pfred", image_tag: str = "tauso/pfred:v1", verbose: bool = True
) -> bool:
    """
    Checks if a specific Docker container is active and running.

    Args:
        container_name: The name of the container to check.
        image_tag: The tag of the image (used for helpful print statements).
        verbose: If True, prints status messages and instructions.

    Returns:
        True if the container is running, False otherwise.
    """
    if verbose:
        print(f"Checking status of container '{container_name}' (Image: {image_tag})...")

    try:
        # Check if the container is running
        result = subprocess.run(
            ["docker", "inspect", "-f", "{{.State.Running}}", container_name],
            capture_output=True,
            text=True,
            check=True,
        )

        if result.stdout.strip() == "true":
            if verbose:
                print(f"✅ [SUCCESS] The container '{container_name}' is active and ready.")
            return True
        else:
            if verbose:
                print(f"⚠️ [WARNING] Container '{container_name}' exists but is STOPPED.")
                print(f"To start it, run: docker start {container_name}")
            return False

    except subprocess.CalledProcessError:
        if verbose:
            print(f"❌ [ERROR] Container '{container_name}' not found.")
            print(f"You need to create it first using your image:")
            print(f"Run: docker run -d -t --name {container_name} {image_tag}")
        return False

    except FileNotFoundError:
        if verbose:
            print("❌ [ERROR] Docker is not installed or not in your system PATH.")
        return False


def score_with_pfred_batch(input_df, container_name="pfred", temp_dir="./pfred_temp_io"):
    """Runs PFRED on a batch of unique sequences via Docker."""
    if input_df.empty:
        return pd.DataFrame()

    os.makedirs(temp_dir, exist_ok=True)
    local_in = os.path.join(temp_dir, "pfred_input.csv")
    local_out = os.path.join(temp_dir, "pfred_results.csv")
    remote_dir = "/home/pfred/scratch"
    remote_in = "pfred_input.csv"

    if os.path.exists(local_out):
        os.remove(local_out)

    # Format for PFRED: [Seq, ID, DummyScore]
    export_df = pd.DataFrame({"col0": input_df["seq"], "col1": input_df["name"], "col2": 1.0})
    export_df.to_csv(local_in, index=False, header=["Seq", "ID", "Score"])

    try:
        subprocess.run(["docker", "exec", container_name, "mkdir", "-p", remote_dir], check=False)
        subprocess.run(["docker", "cp", local_in, f"{container_name}:{remote_dir}/{remote_in}"], check=True)

        cmd = f"""
        cd {remote_dir}
        tr -d '\r' < {remote_in} > {remote_in}.tmp && mv {remote_in}.tmp {remote_in}
        echo -e "15\n21\n100\n1000\n12" > input_params.txt

        export PATH=/home/pfred/bin/R2.6.0/bin:$PATH
        export R_HOME=/home/pfred/bin/R2.6.0/lib64/R
        export PYTHONPATH=/usr/lib64/python2.6/site-packages:$PYTHONPATH

        python2 /home/pfred/scripts/pfred/antisense_predictor.py \
            AOBase \
            /home/pfred/scripts/pfred/AOBase_542seq_cleaned_modelBuilding_Jan2009_15_21_noOutliers.csv \
            c_a_thermo \
            predict \
            ./{remote_in} \
            < input_params.txt > run.log 2>&1

        ls -t *.csv | grep -v "{remote_in}" | head -n 1
        """

        result = subprocess.run(
            ["docker", "exec", container_name, "/bin/bash", "-c", cmd], capture_output=True, text=True
        )

        target_file = result.stdout.strip()
        if not target_file:
            return pd.DataFrame()

        subprocess.run(["docker", "cp", f"{container_name}:{remote_dir}/{target_file}", local_out], check=True)
        subprocess.run(["docker", "exec", container_name, "rm", "-rf", remote_dir], check=False)

        if os.path.exists(local_out):
            return pd.read_csv(local_out)
        return pd.DataFrame()

    except Exception as e:
        print(f"[ERROR] Docker execution failed: {e}")
        return pd.DataFrame()


def populate_pfred(data, container_name="pfred", temp_dir="./pfred_temp_io", chunk_size=5000):
    """
    Main entry point for populating PFRED features.
    Returns the processed dataframe and the list of feature columns added.
    """
    df = data.copy()

    # 1. Standardize Sequences
    df["clean_seq"] = df[SEQUENCE].astype(str).str.upper().str.replace("U", "T")
    df["clean_seq"] = df["clean_seq"].apply(lambda x: re.sub(r"[^ACTG]", "", x))

    # 2. Filter safe sequences (len >= 16)
    safe_mask = df["clean_seq"].str.len() >= 16
    unique_seqs = df.loc[safe_mask, "clean_seq"].unique()

    all_scores = []

    # 3. Batch Process
    for i in range(0, len(unique_seqs), chunk_size):
        batch = unique_seqs[i : i + chunk_size]
        batch_input = pd.DataFrame({"name": [f"s_{k}" for k in range(len(batch))], "seq": batch})

        batch_res = score_with_pfred_batch(batch_input, container_name, temp_dir)

        if not batch_res.empty:
            if "antisense_strand__5_3" in batch_res.columns:
                batch_res = batch_res.rename(columns={"antisense_strand__5_3": "seq"})
            elif "Seq" in batch_res.columns:
                batch_res = batch_res.rename(columns={"Seq": "seq"})

            if "SVMpred" in batch_res.columns:
                all_scores.append(batch_res[["seq", "SVMpred", "PLSpred"]])

    # 4. Merge
    feature_cols = ["PFRED_SVM", "PFRED_PLS"]

    if all_scores:
        score_df = pd.concat(all_scores, ignore_index=True).drop_duplicates(subset=["seq"])
        svm_map = dict(zip(score_df["seq"], score_df["SVMpred"]))
        pls_map = dict(zip(score_df["seq"], score_df["PLSpred"]))

        df["PFRED_SVM"] = df["clean_seq"].map(svm_map)
        df["PFRED_PLS"] = df["clean_seq"].map(pls_map)
    else:
        df["PFRED_SVM"] = None
        df["PFRED_PLS"] = None

    df.drop(columns=["clean_seq"], inplace=True)

    return df, feature_cols
