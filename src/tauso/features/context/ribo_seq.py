import numpy as np
from pathlib import Path

import pandas as pd
import pyBigWig

from ...new_model.consts_dataframe import CANONICAL_GENE

parent = Path(__file__).parent
RIBOSEQ_40S_HUMAN_DATA = parent / 'human_unselected_40S.RiboProElong.bw'

def reduce_values(values, how):
    if values.size == 0:
        return 0.0

    if how == "sum":
        return values.sum()
    elif how == "mean":
        return values.mean()
    elif how == "max":
        return values.max()
    elif how == "nz_mean":
        nz = values[values > 0]
        return nz.mean() if nz.size > 0 else 0.0

def create_ribo_seq_features(organism, aso_df, flanks=(0,10,20,50,100,125, 150), how="mean"):
    if organism != 'human':
        ValueError("Unsupported organism for ribo_seq feature")

    bw = pyBigWig.open(str(RIBOSEQ_40S_HUMAN_DATA))

    if how not in {"sum", "mean", "max", "nz_mean"}:
        raise ValueError(how)

    # 1. Ensure 'strand' is present in required columns
    required = ["chrom", "target_start", "target_end", "gene_start", "gene_end", "index", "strand"]
    df = aso_df.dropna(subset=required).reset_index(drop=True)

    chroms  = df["chrom"].values
    strands = df["strand"].values   # <--- New: Load Strand
    t_start = df["target_start"].astype(int).values
    t_end   = df["target_end"].astype(int).values
    g_start = df["gene_start"].astype(int).values
    g_end   = df["gene_end"].astype(int).values
    indices = df["index"].values

    features = {}

    for i in range(len(df)):
        chrom = chroms[i]
        strand = strands[i]
        chrom_len = bw.chroms()[chrom]

        # Coordinates for this ASO
        ts, te = t_start[i], t_end[i]

        feat = {}

        # ---------------------------------------------------------
        # A. Whole Gene Feature (Global Context)
        # ---------------------------------------------------------
        gs = max(0, g_start[i])
        ge = min(g_end[i], chrom_len)

        if ge > gs:
            values = np.nan_to_num(np.array(bw.values(chrom, gs, ge), dtype=float), nan=0.0)
            feat[f"ribo_gene_{how}"] = reduce_values(values, how)
        else:
            feat[f"ribo_gene_{how}"] = 0.0

        # ---------------------------------------------------------
        # B. Flank Features (Local Context)
        # ---------------------------------------------------------
        for f in flanks:

            # --- 1. The "Old Way" (Combined: Left + Target + Right) ---
            # This window is centered on the ASO and includes the target site itself.
            s_comb = max(0, ts - f)
            e_comb = min(te + f, chrom_len)

            if e_comb > s_comb:
                vals_comb = np.nan_to_num(np.array(bw.values(chrom, s_comb, e_comb), dtype=float), nan=0.0)
                feat[f"ribo_f{f}_{how}"] = reduce_values(vals_comb, how)
            else:
                feat[f"ribo_f{f}_{how}"] = 0.0

            # --- 2. The "New Way" (Split: Upstream vs Downstream) ---
            # These exclude the target site and look strictly at the sides.
            if f > 0:
                # Genomic Left Region [start - f : start]
                s_left = max(0, ts - f)
                e_left = ts
                val_left = 0.0
                if e_left > s_left:
                    arr = np.nan_to_num(np.array(bw.values(chrom, s_left, e_left), dtype=float), nan=0.0)
                    val_left = reduce_values(arr, how)

                # Genomic Right Region [end : end + f]
                s_right = te
                e_right = min(te + f, chrom_len)
                val_right = 0.0
                if e_right > s_right:
                    arr = np.nan_to_num(np.array(bw.values(chrom, s_right, e_right), dtype=float), nan=0.0)
                    val_right = reduce_values(arr, how)

                # Assign based on Strand
                if strand == "+":
                    # + Strand: Left is 5' (Upstream), Right is 3' (Downstream)
                    feat[f"ribo_upstream_{f}_{how}"]   = val_left
                    feat[f"ribo_downstream_{f}_{how}"] = val_right
                else:
                    # - Strand: Right is 5' (Upstream), Left is 3' (Downstream)
                    feat[f"ribo_upstream_{f}_{how}"]   = val_right
                    feat[f"ribo_downstream_{f}_{how}"] = val_left

        features[indices[i]] = feat

    return features


def add_genomic_coordinates(aso_df, mapper):
    out = aso_df.copy()

    chroms, gene_starts, gene_ends, target_starts, target_ends, strands = [], [], [], [], [], []

    for _, row in out.iterrows():
        gene_name = row.get(CANONICAL_GENE)

        # --- CHANGE IS HERE ---
        # We ask for the GENE coordinates, not a specific transcript
        info = mapper.get_gene_coords(gene_name)

        if not info:
            chroms.append(np.nan)
            gene_starts.append(np.nan)
            gene_ends.append(np.nan)
            target_starts.append(np.nan)
            target_ends.append(np.nan)
            strands.append(np.nan)
            continue

        chrom = info["chrom"]
        g_start = info["gene_start"]
        g_end = info["gene_end"]
        strand = info["strand"]

        # ... The math below is unchanged ...
        L = int(row["Location_in_sequence"])
        k = int(row["true_length_of_seq"])
        gene_len = g_end - g_start

        if strand == "+":
            t_start = g_start + (L - 1)
        else:
            t_start = g_start + (gene_len - L - k + 1)

        t_end = t_start + k

        chroms.append(chrom)
        gene_starts.append(g_start)
        gene_ends.append(g_end)
        target_starts.append(t_start)
        target_ends.append(t_end)
        strands.append(strand)
    out["chrom"] = chroms
    out["gene_start"] = gene_starts
    out["gene_end"] = gene_ends
    out["target_start"] = target_starts
    out["target_end"] = target_ends
    out["strand"] = strands

    return out

def populate_ribo_seq(organism, aso_df, flanks=(0, 10, 20, 50, 100, 125, 150), how="mean"):
    """
    Wraps create_ribo_seq_features to return an enriched DataFrame and a list of new column names.

    Args:
        organism (str): Organism name (e.g., 'human').
        aso_df (pd.DataFrame): Input dataframe containing 'index', 'chrom', 'strand', etc.
        flanks (tuple): Flanking distances.
        how (str): Aggregation method ('mean', 'sum', etc.).

    Returns:
        tuple: (pd.DataFrame, list)
            1. The input aso_df with ribo-seq features added as columns.
            2. A list of the names of the newly added features.
    """
    # 1. Generate the features dictionary
    #    Structure: { index_value: {'ribo_feat_1': val, ...}, ... }
    features_dict = create_ribo_seq_features(organism, aso_df, flanks=flanks, how=how)

    # 2. Convert the dictionary to a DataFrame
    #    orient='index' uses the dictionary keys (the ASO indices) as the DataFrame index
    features_df = pd.DataFrame.from_dict(features_dict, orient='index')

    # 3. capture the list of new column names
    new_features_list = list(features_df.columns)

    # 4. Merge back into the original dataframe
    #    We allow for the possibility that 'index' is a column, not the dataframe index.
    if "index" in aso_df.columns:
        # Set the name of the features_df index to 'index' so we can merge on it
        features_df.index.name = "index"

        # Left merge ensures we keep all rows from aso_df, even if feature gen failed for some
        result_df = aso_df.merge(features_df, on="index", how="left")
    else:
        # Fallback: if user provided a DF where the index is already meaningful
        result_df = aso_df.join(features_df, how="left")

    return result_df, new_features_list