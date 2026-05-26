"""Do MFE and accessibility *add* on top of each other?

Picks the best MFE feature and the best accessibility feature (by overall
|Pearson r| vs averaged log-inhibition), then quantifies their joint vs
marginal contribution:

  - feature-feature correlation (how orthogonal the two families are)
  - nested OLS:  y~MFE,  y~Access,  y~MFE+Access   (in-sample + 5-fold CV R^2)
  - delta-R^2 from adding the 2nd feature, and an additivity ratio
    R^2(both) / (R^2(MFE)+R^2(Access))   (~1 => fully additive, <1 => overlap)
  - partial correlations (each feature controlling for the other)
  - per-cohort version (median across cell lines) so it isn't a pooled artifact

Figure: (A) pooled R^2 bars MFE / Access / Both with a perfect-additivity
reference line; (B) per-cohort R^2(both) vs R^2(best single), diagonal = no gain.

Usage:
    python scripts/fold_additivity.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from scripts.make_fold_graphs import (  # noqa: E402
    ACCESS_FEATURES,
    FEAT_DIR,
    MERGED,
    MFE_FEATURES,
    OUT,
    TARGET,
    load_feature,
)


def ols_r2(X, y):
    """In-sample R^2 of OLS y ~ X (X is (n,k)), with intercept."""
    X1 = np.column_stack([np.ones(len(X)), X])
    beta, *_ = np.linalg.lstsq(X1, y, rcond=None)
    yhat = X1 @ beta
    ss_res = float(((y - yhat) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    return 1.0 - ss_res / ss_tot


def cv_r2(X, y, k=5, seed=0):
    """k-fold CV R^2 (sklearn r2_score convention: test-set mean baseline)."""
    n = len(y)
    idx = np.arange(n)
    np.random.default_rng(seed).shuffle(idx)
    folds = np.array_split(idx, k)
    out = []
    for i in range(k):
        te = folds[i]
        tr = np.concatenate([folds[j] for j in range(k) if j != i])
        Xtr = np.column_stack([np.ones(len(tr)), X[tr]])
        beta, *_ = np.linalg.lstsq(Xtr, y[tr], rcond=None)
        Xte = np.column_stack([np.ones(len(te)), X[te]])
        yhat = Xte @ beta
        ss_res = float(((y[te] - yhat) ** 2).sum())
        ss_tot = float(((y[te] - y[te].mean()) ** 2).sum())
        out.append(1.0 - ss_res / ss_tot)
    return float(np.mean(out))


def residualize(a, b):
    """Residual of a after regressing on b (with intercept)."""
    B = np.column_stack([np.ones(len(b)), b])
    beta, *_ = np.linalg.lstsq(B, a, rcond=None)
    return a - B @ beta


def pick_best(feature_names, y_series):
    """Return (name, |pearson r|) of the feature most correlated with y_series."""
    best, best_abs = None, -1.0
    for name in feature_names:
        fs = load_feature(name)
        if fs is None:
            continue
        df = pd.concat([fs, y_series], axis=1).dropna()
        if len(df) < 100 or df[name].std() < 1e-12:
            continue
        r = abs(pearsonr(df[name].to_numpy(float), df.iloc[:, 1].to_numpy(float))[0])
        if r > best_abs:
            best, best_abs = name, r
    return best, best_abs


def main():
    target = pd.read_csv(TARGET).rename(columns={"Inhibition(%)": "inhibition_percent"})
    target = target.set_index("index_oligo")["inhibition_percent"]
    target = target[target > 0]
    logy = np.log(target)
    logy.name = "logy"

    mfe_name, mfe_r = pick_best(MFE_FEATURES, logy)
    acc_name, acc_r = pick_best(ACCESS_FEATURES, logy)
    print(f"best MFE   : {mfe_name}  (|r|={mfe_r:.3f})")
    print(f"best ACCESS: {acc_name}  (|r|={acc_r:.3f})")

    mfe = load_feature(mfe_name)
    acc = load_feature(acc_name)
    df = pd.concat([mfe, acc, logy], axis=1).dropna()
    m = df[mfe_name].to_numpy(float)
    a = df[acc_name].to_numpy(float)
    y = df["logy"].to_numpy(float)
    n = len(df)
    print(f"\naligned N = {n:,}")

    # how independent are the two features?
    ff_p = pearsonr(m, a)[0]
    ff_s = spearmanr(m, a)[0]
    print(f"\nfeature-feature corr:  Pearson {ff_p:+.3f}   Spearman {ff_s:+.3f}")

    # nested models
    r2_m = ols_r2(m[:, None], y)
    r2_a = ols_r2(a[:, None], y)
    r2_both = ols_r2(np.column_stack([m, a]), y)
    cv_m = cv_r2(m[:, None], y)
    cv_a = cv_r2(a[:, None], y)
    cv_both = cv_r2(np.column_stack([m, a]), y)

    print("\n--- pooled (averaged log-inhibition) ---")
    print(f"{'model':<18}{'R2 (in-sample)':>16}{'R2 (5-fold CV)':>16}")
    print(f"{'MFE only':<18}{r2_m:>16.4f}{cv_m:>16.4f}")
    print(f"{'Access only':<18}{r2_a:>16.4f}{cv_a:>16.4f}")
    print(f"{'MFE + Access':<18}{r2_both:>16.4f}{cv_both:>16.4f}")
    print(f"\ndelta-R2 add Access to MFE : {r2_both - r2_m:+.4f}")
    print(f"delta-R2 add MFE to Access : {r2_both - r2_a:+.4f}")
    print(f"additivity ratio R2(both)/(R2_MFE+R2_Access) = {r2_both/(r2_m+r2_a):.3f}  (1.0=fully additive)")

    # partial correlations
    pc_access = pearsonr(residualize(a, m), residualize(y, m))[0]
    pc_mfe = pearsonr(residualize(m, a), residualize(y, a))[0]
    print(f"\npartial corr Access | MFE  : {pc_access:+.3f}  (vs raw {pearsonr(a,y)[0]:+.3f})")
    print(f"partial corr MFE | Access  : {pc_mfe:+.3f}  (vs raw {pearsonr(m,y)[0]:+.3f})")

    # per-cohort
    merged = pd.read_parquet(MERGED)[["Cell_line", "inhibition_percent"]]
    joined = merged.join(mfe, how="inner").join(acc, how="inner")
    joined = joined[joined["inhibition_percent"] > 0]
    joined["logy"] = np.log(joined["inhibition_percent"])
    rows = []
    for cell, sub in joined.groupby("Cell_line"):
        s = sub[[mfe_name, acc_name, "logy"]].dropna()
        if len(s) < 400 or s[mfe_name].std() < 1e-12 or s[acc_name].std() < 1e-12:
            continue
        ym = s["logy"].to_numpy(float)
        mm = s[mfe_name].to_numpy(float)
        aa = s[acc_name].to_numpy(float)
        rows.append((
            cell, len(s),
            ols_r2(mm[:, None], ym),
            ols_r2(aa[:, None], ym),
            ols_r2(np.column_stack([mm, aa]), ym),
        ))
    coh = pd.DataFrame(rows, columns=["cell", "n", "r2_mfe", "r2_acc", "r2_both"])
    coh["r2_best_single"] = coh[["r2_mfe", "r2_acc"]].max(axis=1)
    print(f"\n--- per-cohort ({len(coh)} cell lines, n>=400) ---")
    print(f"median R2  MFE={coh.r2_mfe.median():.4f}  Access={coh.r2_acc.median():.4f}  "
          f"Both={coh.r2_both.median():.4f}")
    print(f"median delta-R2 Both - best single = {(coh.r2_both - coh.r2_best_single).median():+.4f}")
    print(f"cohorts where Both > best single   = {(coh.r2_both > coh.r2_best_single + 1e-6).sum()}/{len(coh)}")
    coh.sort_values("r2_both", ascending=False).to_csv(OUT / "fold_additivity_per_cohort.csv", index=False)

    # ---- figure ----
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))

    labels = ["MFE\nonly", "Access\nonly", "MFE +\nAccess"]
    vals = [r2_m, r2_a, r2_both]
    cvs = [cv_m, cv_a, cv_both]
    xpos = np.arange(3)
    ax1.bar(xpos - 0.18, vals, width=0.36, color=["#2ca02c", "#d62728", "#1f77b4"], label="in-sample")
    ax1.bar(xpos + 0.18, cvs, width=0.36, color=["#2ca02c", "#d62728", "#1f77b4"], alpha=0.45, label="5-fold CV")
    ax1.axhline(r2_m + r2_a, ls="--", color="k", lw=1.2, label="perfect additivity (R²_MFE+R²_Access)")
    for xi, v in zip(xpos, vals):
        ax1.text(xi - 0.18, v + 0.0004, f"{v:.4f}", ha="center", va="bottom", fontsize=8)
    ax1.set_xticks(xpos)
    ax1.set_xticklabels(labels)
    ax1.set_ylabel("R²  (averaged log-inhibition)")
    ax1.set_title(f"(A) pooled  N={n:,}\nfeat-feat r={ff_p:+.2f}  →  ΔR²(+Access)={r2_both-r2_m:+.4f}")
    ax1.legend(fontsize=8, loc="upper left")

    ax2.scatter(coh["r2_best_single"], coh["r2_both"], s=22, alpha=0.7, color="#1f77b4")
    lim = max(coh["r2_both"].max(), coh["r2_best_single"].max()) * 1.1
    ax2.plot([0, lim], [0, lim], ls="--", color="k", lw=1, label="no gain (y=x)")
    ax2.set_xlim(0, lim)
    ax2.set_ylim(0, lim)
    ax2.set_xlabel("R² best single feature (per cohort)")
    ax2.set_ylabel("R² MFE + Access (per cohort)")
    ax2.set_title(f"(B) per-cohort: {(coh.r2_both > coh.r2_best_single+1e-6).sum()}/{len(coh)} above diagonal")
    ax2.legend(fontsize=8, loc="upper left")

    fig.suptitle(f"Do MFE & accessibility add?   {mfe_name}  +  {acc_name}", fontsize=12)
    fig.tight_layout()
    out_path = OUT / "fold_mfe_access_additivity.png"
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()
