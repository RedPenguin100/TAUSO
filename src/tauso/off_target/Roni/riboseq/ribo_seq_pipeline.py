import os
import pandas as pd
import numpy as np
from pyprojroot import here
import pyBigWig

from notebooks.consts import NOTEBOOK_PATH
from src.tauso.off_target.Roni.off_target_pipeline.off_target_functions import parse_gtf
from notebooks.consts import SAVED_FEATURES

PROJECT_ROOT = here()
DATA_DIR = os.path.join(PROJECT_ROOT, "src", "tauso", "off_target", "Roni", "data")

gtf_path = os.path.join(DATA_DIR, "gencode.v34.chr_patch_hapl_scaff.annotation.gtf")
# Genomic annotation data
annotations = parse_gtf(gtf_path, "../outputs/gtf_pickle_output_file_path.pkl")
# ASO data
aso_df = pd.read_csv(NOTEBOOK_PATH / "data/data_asoptimizer_updated.csv")
# Ribo-seq data
bw = pyBigWig.open("human_unselected_40S.RiboProElong.bw")

def build_gene_lookup_from_transcripts(annotations):
    """
    Build:
    gene_name -> {chrom, gene_start, gene_end, strand}

    Resolves strand conflicts by choosing the longest transcript.
    """

    gene_to_transcripts = {}

    for tx_id, info in annotations.items():
        gene = info["gene_name"]

        start0 = info["start"] - 1  # 0-based
        end0   = info["end"]
        length = end0 - start0

        record = {
            "chrom": info["chrom"],
            "gene_start": start0,
            "gene_end": end0,
            "strand": info["strand"],
            "length": length
        }

        gene_to_transcripts.setdefault(gene, []).append(record)

    gene_lookup = {}

    for gene, records in gene_to_transcripts.items():
        # choose longest transcript
        longest = max(records, key=lambda r: r["length"])

        gene_lookup[gene] = {
            "chrom": longest["chrom"],
            "gene_start": longest["gene_start"],
            "gene_end": longest["gene_end"],
            "strand": longest["strand"]
        }

    return gene_lookup

def add_genomic_coordinates(aso_df, annotations):
    """
    Adds genomic coordinates to aso_df:
    chrom, gene_start, gene_end, target_start, target_end

    Rows with unmappable genes get NaN values.
    """

    gene_lookup = build_gene_lookup_from_transcripts(annotations)

    out = aso_df.copy()

    chroms = []
    gene_starts = []
    gene_ends = []
    target_starts = []
    target_ends = []

    for _, row in out.iterrows():
        gene = row.get("Canonical Gene Name")

        # Handle controls / missing genes
        if not isinstance(gene, str) or gene not in gene_lookup:
            chroms.append(np.nan)
            gene_starts.append(np.nan)
            gene_ends.append(np.nan)
            target_starts.append(np.nan)
            target_ends.append(np.nan)
            continue

        info = gene_lookup[gene]
        chrom = info["chrom"]
        g_start = info["gene_start"]
        g_end = info["gene_end"]
        strand = info["strand"]

        L = int(row["Location_in_sequence"])   # 1-based
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

    out["chrom"] = chroms
    out["gene_start"] = gene_starts
    out["gene_end"] = gene_ends
    out["target_start"] = target_starts
    out["target_end"] = target_ends

    return out

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

def create_ribo_seq_features(bw, aso_df, flanks=(0,10,20,50,100), how="mean"):

    if how not in {"sum", "mean", "max", "nz_mean"}:
        raise ValueError(how)

    required = ["chrom", "target_start", "target_end", "gene_start", "gene_end", "index"]
    df = aso_df.dropna(subset=required).reset_index(drop=True)

    chroms  = df["chrom"].values
    t_start = df["target_start"].astype(int).values
    t_end   = df["target_end"].astype(int).values
    g_start = df["gene_start"].astype(int).values
    g_end   = df["gene_end"].astype(int).values
    indices = df["index"].values

    features = {}

    for i in range(len(df)):
        chrom = chroms[i]
        chrom_len = bw.chroms()[chrom]

        feat = {}

        # ---------- gene ----------
        gs = max(0, g_start[i])
        ge = min(g_end[i], chrom_len)

        if ge > gs:
            values = np.nan_to_num(
                np.array(bw.values(chrom, gs, ge), dtype=float),
                nan=0.0
            )
            feat[f"ribo_gene_{how}"] = reduce_values(values, how)
        else:
            feat[f"ribo_gene_{how}"] = 0.0

        # ---------- flanks ----------
        for f in flanks:
            s = max(0, t_start[i] - f)
            e = min(t_end[i] + f, chrom_len)

            if e > s:
                values = np.nan_to_num(
                    np.array(bw.values(chrom, s, e), dtype=float),
                    nan=0.0
                )
                val = reduce_values(values, how)
            else:
                val = 0.0

            feat[f"ribo_f{f}_{how}"] = val

        features[indices[i]] = feat

    return features

def save_ribo_seq_features(feature_dict, feature_dir):
    """
    Save ribo-seq features to CSV files.
    One CSV per feature:
        index, feature_value
    """
    os.makedirs(feature_dir, exist_ok=True)

    # Collect all feature names
    feature_names = set()
    for feats in feature_dict.values():
        feature_names.update(feats.keys())

    # Build one dataframe per feature
    for feature in sorted(feature_names):
        rows = []

        for aso_idx, feats in feature_dict.items():
            rows.append({
                "index": aso_idx,
                feature: feats.get(feature, np.nan)  # <-- FIX
            })

        df = pd.DataFrame(rows)
        df.sort_values("index", inplace=True)

        out_path = os.path.join(feature_dir, f"{feature}.csv")
        df.to_csv(out_path, index=False)

def generate_riboseq_features(bw=bw,
                              aso_df=aso_df,
                              annotations=annotations,
                              flanks=(0,10,20,50,100),
                              how="mean"
                              ):
    new = add_genomic_coordinates(aso_df, annotations)
    features = create_ribo_seq_features(bw, new, flanks, how)
    save_ribo_seq_features(features, feature_dir=SAVED_FEATURES)