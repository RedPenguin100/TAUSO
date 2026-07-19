"""Shared loader for the model-comparison figures.

TAUSO vs competitor scores on the held-out test split (strict gapmers, rows OligoAI also scored).
TAUSO is the shipped champion, trained on train+val -- it has not seen the test rows -- and scored
once on the test split. Everything is keyed by index_oligo: competitor scores come from the
feature-store competition shards, TAUSO's test predictions from tauso_test_pred.parquet.
"""
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[3]          # notebooks/graphs/models -> repo root
sys.path.insert(0, str(REPO)); sys.path.insert(0, str(REPO / "src"))
from notebooks.data.OligoAI.parse_chemistry import assign_chemistry    # noqa: E402
from tauso.data.consts import CHEMICAL_PATTERN                         # noqa: E402
from tauso.data.data import get_data_dir                               # noqa: E402

PATCHES = str(Path(get_data_dir()) / "features" / "oligo" / "_competitors_v15")
AVG = str(REPO / "notebooks/data/OligoAI/raw_data/aso_inhibitions_with_canonical_gene_processed_averaged.csv.gz")
IDX_CSV = str(REPO / "notebooks/data/OligoAI/raw_data/aso_inhibitions_with_canonical_gene_indexed.csv.gz")
TAUSO_PRED = str(Path(__file__).resolve().parent / "tauso_test_pred.parquet")

IDX, INH, CID = "index_oligo", "inhibition_percent", "custom_id"
GENE, CELL = "canonical_gene_name", "cell_line"
STRICT = {"MMMMMddddddddddMMMMM", "CCCddddddddddCCC"}   # 5-10-5 MOE / 3-10-3 cEt gapmers
COMP = {"OligoAI": "oligo_ai_score", "OligoWalk·intra": "OW_Intra_Oligo", "OligoWalk·Tm": "OW_Tm",
        "OligoWalk·duplex": "OW_Duplex", "OligoWalk": "OW_Overall", "miRanda·E": "miranda_energy",
        "miRanda": "miranda_score", "sfold": "sfold_accessibility", "PFRED": "PFRED_PLS", "PFRED·SVM": "PFRED_SVM"}


def _shard(stem):
    p = f"{PATCHES}/{stem}.parquet"
    d = pd.read_parquet(p) if os.path.exists(p) else pd.read_csv(f"{PATCHES}/{stem}.csv")
    return d[[IDX, stem]]


def load():
    avg = pd.read_csv(AVG, usecols=[IDX, INH, CID, "split", CHEMICAL_PATTERN, GENE, CELL])
    avg["is_strict"] = avg[CHEMICAL_PATTERN].isin(STRICT)
    df = avg[[IDX, INH, CID, GENE, CELL, "split", "is_strict"]]
    present = {}
    for lbl, stem in COMP.items():
        try:
            df = df.merge(_shard(stem), on=IDX, how="left"); present[lbl] = stem
        except FileNotFoundError:
            pass
    df = df.merge(pd.read_parquet(TAUSO_PRED), on=IDX, how="left"); present["TAUSO"] = "tauso"
    df = df[(df["split"] == "test") & df["is_strict"] & df[INH].notna()
            & df["oligo_ai_score"].notna()].reset_index(drop=True)   # apples-to-apples vs OligoAI
    df["gxc"] = df[GENE].astype(str) + "|" + df[CELL].astype(str)
    return df, present


def med(df, score, grp):
    """Median per-cohort Spearman of `score` vs inhibition (cohorts with >1 unique score, >=3 rows)."""
    d = df[[score, INH, grp]].dropna()
    rs = [spearmanr(g[score], g[INH]).correlation for _, g in d.groupby(grp)
          if g[score].nunique() > 1 and len(g) >= 3]
    rs = [r for r in rs if np.isfinite(r)]
    return float(np.median(rs)) if rs else np.nan


def orient(df, present):
    """Sign per method: flip a score stored 'backwards' (anti-correlated in BOTH groupings). TAUSO stays +1."""
    sign = {}
    for lbl, col in present.items():
        c, g = med(df, col, CID), med(df, col, "gxc")
        sign[lbl] = -1.0 if (np.isfinite(c) and np.isfinite(g) and c < 0 and g < 0) else 1.0
    return sign
