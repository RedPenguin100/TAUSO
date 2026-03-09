from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig

from ...new_model.consts_dataframe import CANONICAL_GENE, SENSE_LENGTH, SENSE_START

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


def calculate_ribo_seq_row(row, bw, flanks, how):
    chrom = row["chrom"]
    strand = row["strand"]
    chrom_len = bw.chroms()[chrom]

    ts, te = int(row["target_start"]), int(row["target_end"])
    gs, ge = int(row["gene_start"]), int(row["gene_end"])

    feat = {}

    # Global Context
    gs = max(0, gs)
    ge = min(ge, chrom_len)

    if ge > gs:
        values = np.nan_to_num(np.array(bw.values(chrom, gs, ge), dtype=float), nan=0.0)
        feat[f"ribo_gene_{how}"] = reduce_values(values, how)
    else:
        feat[f"ribo_gene_{how}"] = 0.0

    # Local Context
    for f in flanks:
        s_comb = max(0, ts - f)
        e_comb = min(te + f, chrom_len)

        if e_comb > s_comb:
            vals_comb = np.nan_to_num(np.array(bw.values(chrom, s_comb, e_comb), dtype=float), nan=0.0)
            feat[f"ribo_f{f}_{how}"] = reduce_values(vals_comb, how)
        else:
            feat[f"ribo_f{f}_{how}"] = 0.0

        if f > 0:
            s_left = max(0, ts - f)
            e_left = ts
            val_left = 0.0
            if e_left > s_left:
                arr = np.nan_to_num(np.array(bw.values(chrom, s_left, e_left), dtype=float), nan=0.0)
                val_left = reduce_values(arr, how)

            s_right = te
            e_right = min(te + f, chrom_len)
            val_right = 0.0
            if e_right > s_right:
                arr = np.nan_to_num(np.array(bw.values(chrom, s_right, e_right), dtype=float), nan=0.0)
                val_right = reduce_values(arr, how)

            if strand == "+":
                feat[f"ribo_upstream_{f}_{how}"] = val_left
                feat[f"ribo_downstream_{f}_{how}"] = val_right
            else:
                feat[f"ribo_upstream_{f}_{how}"] = val_right
                feat[f"ribo_downstream_{f}_{how}"] = val_left

    return feat


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
        L = int(row[SENSE_START])
        k = int(row[SENSE_LENGTH])
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
