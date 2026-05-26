"""Generate article-style graphs for the (re-parametrized) fold features.

For each feature:
  - 2-panel figure (A: scatter vs log-inhibition + regression line;
    B: binned trend, mean log-inhibition per quantile bin, SEM bars + 95% CI),
    matching the Figure X1/X2 style.
Plus one summary bar chart: median Spearman across cell-line cohorts per feature.

Inhibition for the scatter uses the averaged target table (one row per oligo,
~N=30k, as in the paper). Per-cohort Spearman uses merged_curated (per cell line).

Usage:
    python scripts/make_fold_graphs.py --family mfe
    python scripts/make_fold_graphs.py --family access
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import linregress, pearsonr, spearmanr

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
FEAT_DIR = REPO / "notebooks" / "saved_features_oligo"
TARGET = REPO / "notebooks" / "data" / "OligoAI" / "raw_data" / "aso_inhibitions_with_canonical_gene_processed_averaged.csv.gz"
MERGED = REPO / "notebooks" / "models" / "FeatureAnalysis" / "cache" / "merged_curated.parquet"
OUT = REPO / "notebooks" / "graphs" / "features"
OUT.mkdir(parents=True, exist_ok=True)

# current-default fold features (must match populate_fold.py)
MFE_FEATURES = [
    "mfe_win45_flank30_step7", "mfe_win50_flank30_step7", "mfe_win55_flank30_step7",
    "mfe_win60_flank30_step7", "mfe_win65_flank30_step7", "mfe_win55_flank30_step15",
    "mfe_win65_flank30_step15", "mfe_win55_flank60_step7", "mfe_win55_flank120_step7",
    "mfe_win65_flank120_step7", "mfe_win20_flank30_step4", "mfe_win30_flank30_step7",
    "mfe_win90_flank120_step15", "mfe_win45_flank240_step7", "mfe_win90_flank240_step15",
]
ACCESS_FEATURES = [
    "access_120flank_6access_4-6seed_sizes", "access_120flank_8access_4-6-8seed_sizes",
    "access_120flank_10access_4-6-8seed_sizes", "access_120flank_13access_4-6-8seed_sizes",
    "access_120flank_13access_13-26-39seed_sizes", "access_120flank_20access_4-6-8seed_sizes",
    "access_120flank_39access_13-26-39seed_sizes",
]


def load_feature(name):
    p = FEAT_DIR / f"{name}.csv"
    if not p.exists():
        return None
    s = pd.read_csv(p).set_index("index_oligo").iloc[:, 0]
    s.name = name
    return s


def two_panel(feat_name, x, y, out_path):
    """x=feature values, y=log inhibition (aligned, NaN-free)."""
    rho, rho_p = spearmanr(x, y)
    r, r_p = pearsonr(x, y)
    lr = linregress(x, y)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    # Panel A: scatter + regression
    ax1.scatter(x, y, s=4, alpha=0.08, color="#3b6", rasterized=True)
    xs = np.linspace(np.nanmin(x), np.nanmax(x), 100)
    ax1.plot(xs, lr.intercept + lr.slope * xs, color="black", lw=2)
    ax1.set_xlabel(feat_name)
    ax1.set_ylabel("log(inhibition)")
    ax1.set_title(f"(A) raw + linear fit  (N={len(x):,})")

    # Panel B: binned trend with SEM + 95% CI
    nbins = 20
    edges = np.quantile(x, np.linspace(0, 1, nbins + 1))
    edges = np.unique(edges)
    idx = np.clip(np.digitize(x, edges[1:-1]), 0, len(edges) - 2)
    centers, means, sems = [], [], []
    for b in range(len(edges) - 1):
        m = idx == b
        if m.sum() < 5:
            continue
        centers.append(x[m].mean())
        means.append(y[m].mean())
        sems.append(y[m].std(ddof=1) / np.sqrt(m.sum()))
    centers, means, sems = map(np.array, (centers, means, sems))
    ax2.errorbar(centers, means, yerr=sems, fmt="o", color="#1f77b4", capsize=3, label="mean ± SEM")
    ax2.fill_between(centers, means - 1.96 * sems, means + 1.96 * sems, alpha=0.2, color="#1f77b4", label="95% CI")
    ax2.set_xlabel(feat_name)
    ax2.set_ylabel("mean log(inhibition)")
    ax2.set_title("(B) binned trend")
    ax2.legend()

    fig.suptitle(f"{feat_name}   Spearman ρ={rho:.3f} (p={rho_p:.1e})   Pearson r={r:.3f}",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return dict(feature=feat_name, n=len(x), pearson_r=r, spearman_rho=rho, spearman_p=rho_p)


def cohort_median_spearman(merged, feat_series, feat_name, min_cohort=400):
    df = merged.join(feat_series, how="inner")
    rhos = []
    for cell, sub in df.groupby("Cell_line"):
        s = sub[[feat_name, "inhibition_percent"]].dropna()
        if len(s) < min_cohort or s[feat_name].std() < 1e-12:
            continue
        rhos.append(spearmanr(s[feat_name], s["inhibition_percent"]).statistic)
    rhos = np.array(rhos)
    if len(rhos) == 0:
        return np.nan, 0
    return float(np.median(rhos)), len(rhos)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", choices=["mfe", "access"], required=True)
    args = ap.parse_args()
    feats = MFE_FEATURES if args.family == "mfe" else ACCESS_FEATURES

    target = pd.read_csv(TARGET).rename(columns={"Inhibition(%)": "inhibition_percent"})
    target = target.set_index("index_oligo")["inhibition_percent"]
    merged = pd.read_parquet(MERGED)[["Cell_line", "inhibition_percent"]]

    stats = []
    for name in feats:
        fs = load_feature(name)
        if fs is None:
            print(f"  SKIP (missing): {name}")
            continue
        # scatter vs averaged log-inhibition
        df = pd.concat([fs, target], axis=1).dropna()
        df = df[df["inhibition_percent"] > 0]
        if len(df) < 100 or df[name].std() < 1e-12:
            print(f"  SKIP (degenerate/all-NaN, {len(df)} valid rows): {name}")
            stats.append(dict(feature=name, n=len(df), pearson_r=np.nan,
                              spearman_rho=np.nan, spearman_p=np.nan,
                              cohort_median_spearman=np.nan, n_cohorts=0))
            continue
        x = df[name].to_numpy(float)
        y = np.log(df["inhibition_percent"].to_numpy(float))
        st = two_panel(name, x, y, OUT / f"{name}__vs_loginhibition.png")
        med_rho, ncoh = cohort_median_spearman(merged, fs, name)
        st["cohort_median_spearman"] = med_rho
        st["n_cohorts"] = ncoh
        stats.append(st)
        print(f"  {name}: spearman={st['spearman_rho']:+.3f}  cohort_median={med_rho:+.3f} ({ncoh} cohorts)")

    sdf = pd.DataFrame(stats).sort_values("cohort_median_spearman", key=lambda s: s.abs(), ascending=False)
    sdf.to_csv(OUT / f"{args.family}_fold_stats.csv", index=False)

    # summary bar chart
    fig, ax = plt.subplots(figsize=(10, max(4, 0.45 * len(sdf))))
    ax.barh(sdf["feature"], sdf["cohort_median_spearman"], color="#d62728" if args.family == "access" else "#2ca02c")
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("median Spearman across cohorts")
    ax.set_title(f"{args.family.upper()} fold features — per-cohort median Spearman (new parametrization)")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(OUT / f"{args.family}_cohort_median_spearman.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    print(f"\nWrote {len(stats)} per-feature figures + summary to {OUT}")


if __name__ == "__main__":
    main()
