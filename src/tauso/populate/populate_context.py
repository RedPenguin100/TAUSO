from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import pyBigWig

from ..data.consts import CANONICAL_GENE, CELL_LINE_DEPMAP
from ..features.context.ribo_seq import RIBOSEQ_40S_HUMAN_DATA, calculate_ribo_seq_row


def populate_ribo_seq(organism, aso_df, flanks=(0, 10, 20, 50, 100, 125, 150), how="mean"):
    if organism != "human":
        raise ValueError("Unsupported organism for ribo_seq feature")

    if how not in {"sum", "mean", "max", "nz_mean"}:
        raise ValueError(how)

    bw = pyBigWig.open(str(RIBOSEQ_40S_HUMAN_DATA))

    # Pre-initialize columns with NaN to ensure structure exists even if rows fail
    # We define the column names first based on flanks and 'how'
    new_features_list = [f"ribo_gene_{how}"]
    for f in flanks:
        new_features_list.append(f"ribo_f{f}_{how}")
        if f > 0:
            new_features_list.append(f"ribo_upstream_{f}_{how}")
            new_features_list.append(f"ribo_downstream_{f}_{how}")

    # Initialize new columns in the original dataframe
    for col in new_features_list:
        aso_df[col] = np.nan

    # Iterate by index to assign directly in place
    for idx, row in aso_df.iterrows():
        try:
            feat_dict = calculate_ribo_seq_row(row, bw, flanks, how)

            # Assign values directly to the specific index
            for key, val in feat_dict.items():
                aso_df.at[idx, key] = val

        except Exception as e:
            print(f"Row {idx} failed: {e}")
            continue

    return aso_df, new_features_list


def populate_mrna_expression(
    df: pd.DataFrame, expression_dict: Dict[str, pd.DataFrame]
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Enriches the main dataframe with mRNA expression data for the target gene
    and the specific RNASEH1 gene.
    """

    # --- MINIMAL FIX 1: Warn and drop to prevent column duplication ---
    existing_cols = [c for c in ["target_expression", "rnase_expression"] if c in df.columns]
    if existing_cols:
        print(f"WARNING: Dropping existing columns to prevent alignment breaks: {existing_cols}")
        df = df.drop(columns=existing_cols)
    # ------------------------------------------------------------------

    # 1. Flatten the dictionary into one Master DataFrame
    # ---------------------------------------------------------
    dfs_to_concat = []

    # We iterate through the dictionary to add the cell line key as a column
    for depmap_id, t_df in expression_dict.items():
        # Defensive copy to avoid SettingWithCopy warnings on the source DFs
        temp = t_df[["Gene", "expression_norm"]].copy()
        temp[CELL_LINE_DEPMAP] = depmap_id
        dfs_to_concat.append(temp)

    expression_master = pd.concat(dfs_to_concat, ignore_index=True)

    # --- MINIMAL FIX 2: Prevent row multiplication during merge ---
    expression_master = expression_master.drop_duplicates(subset=[CELL_LINE_DEPMAP, "Gene"])
    # --------------------------------------------------------------

    # 2. Get 'target_expression' (General Merge)
    # ---------------------------------------------------------
    enhanced_df = df.merge(
        expression_master,
        left_on=[CELL_LINE_DEPMAP, CANONICAL_GENE],
        right_on=[CELL_LINE_DEPMAP, "Gene"],
        how="left",
    )

    # Rename and clean up
    enhanced_df = enhanced_df.rename(columns={"expression_norm": "target_expression"})
    enhanced_df = enhanced_df.drop(columns=["Gene"])

    # 3. Get 'rnase_expression' (Specific 'RNASEH1' Merge)
    # ---------------------------------------------------------
    # Filter master list for RNASEH1
    rnase_vals = expression_master[expression_master["Gene"] == "RNASEH1"]

    # Merge only on cell line ID
    enhanced_df = enhanced_df.merge(
        rnase_vals[[CELL_LINE_DEPMAP, "expression_norm"]],
        on=CELL_LINE_DEPMAP,
        how="left",
        suffixes=("", "_rnase"),  # Safety against collisions
    )

    # Rename
    enhanced_df = enhanced_df.rename(columns={"expression_norm": "rnase_expression"})

    # Define the list of features added
    new_features = ["target_expression", "rnase_expression"]

    return enhanced_df, new_features


def populate_transfection(data):
    """
    One-hot encodes the specified column into 0s and 1s, ensures all expected
    features exist (even if absent from the data), and merges them.
    """
    # The master list of expected categories
    features = ["Electroporation", "Gymnosis", "Lipofection", "Other"]

    # 1. Create binary columns as 0/1 for whatever actually exists in the data
    binary_features = pd.get_dummies(data["transfection_method"], dtype=int)

    # 2. Force the dataframe to have exactly our expected columns, filling missing ones with 0
    binary_features = binary_features.reindex(columns=features, fill_value=0)

    # 3. Join them back to your original dataframe
    data = pd.concat([data, binary_features], axis=1)

    return data, features
