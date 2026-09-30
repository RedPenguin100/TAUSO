"""The parsimony curve on its own scale: what removing features costs, CV above and TEST below.

Companion to plot_vs_oligoai.py, which shares a y-axis with the competitor line and so hides
our own drift. Error bars are the standard error of the mean (CV sd/sqrt(15), TEST sd/sqrt(3));
the raw fold-to-fold sd is an order of magnitude larger and describes how much folds differ
from each other, not how precisely each rung is measured.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ===== TWEAK ME =====
L = Path(__file__).resolve().parent
OUT = L / "parsimony_curve.png"
C_CV, C_TE = "#1f5f8b", "#b4491f"
FIGSIZE = (16.0, 6.4)
TICKS = (679, 500, 400, 300, 200, 150, 100, 60, 40, 25, 15, 10)
N_CV_FITS, N_TEST_SEEDS = 15, 3
# =====================

other = pd.read_csv(L / "other_models_test.csv")
m = pd.read_csv(L / "curve.csv").merge(pd.read_csv(L / "test_curve.csv"),
                                       on=["step", "n"], suffixes=("", "_t"))
metrics = [("rmse", "RMSE (centered target)", True),
           ("exp_med", "exp_med (within custom_id)", False),
           ("exp_mean", "exp_mean (within custom_id)", False),
           ("gxc_med", "gxc_med (within gene x cell)", False)]

fig, axes = plt.subplots(2, len(metrics), figsize=FIGSIZE, sharex=True)
for col, (key, label, lower_better) in enumerate(metrics):
    for row, (prefix, color, nrep, tag) in enumerate(
            [("cv", C_CV, N_CV_FITS, f"CV ({N_CV_FITS} fits)"),
             ("test", C_TE, N_TEST_SEEDS, f"TEST ({N_TEST_SEEDS} seeds)")]):
        ax = axes[row][col]
        y = m[f"{prefix}_{key}"]
        sem = m[f"{prefix}_{key}_sd"] / np.sqrt(nrep)
        ax.errorbar(m.n, y, yerr=sem, fmt="o-", ms=3.2, lw=0.9, elinewidth=0.7,
                    capsize=1.5, color=color)
        if prefix == "test" and key == "exp_med":
            for o in other.itertuples():
                ax.axhline(o.test_exp_med, color=o.color, ls=":", lw=1.3)
        pick = y.idxmin() if lower_better else y.idxmax()
        ax.plot(m.n[pick], y[pick], marker="*", ms=14, color=color, lw=0, zorder=5)
        ax.set_xscale("log")
        ax.set_xlim(m.n.max() * 1.07, m.n.min() * 0.93)
        ticks = [t for t in TICKS if m.n.min() * 0.95 <= t <= m.n.max() * 1.05]
        ax.set_xticks(ticks); ax.set_xticklabels([str(t) for t in ticks]); ax.minorticks_off()
        ax.grid(alpha=0.22, lw=0.5)
        ax.set_ylabel(label if col == 0 else "", fontsize=8)
        ax.set_title(f"{tag}  |  {label}\nbest n={int(m.n[pick])} ({y[pick]:.4f})", fontsize=8.5)
        if row == 1:
            ax.set_xlabel("features in the model (descending)", fontsize=9)

fig.suptitle(f"TAUSO parsimony descent, {len(m)} rungs: n={int(m.n.max())} -> {int(m.n.min())}"
             f"   (error bars = standard error of the mean)", fontsize=11)
fig.tight_layout()
fig.savefig(OUT, dpi=200)
print(f"wrote {OUT}")
for key, _, lower in metrics:
    for p, nrep in (("cv", N_CV_FITS), ("test", N_TEST_SEEDS)):
        c = f"{p}_{key}"
        print(f"{c:16s} {m[c].iloc[0]:.4f} -> {m[c].iloc[-1]:.4f}  "
              f"(delta {m[c].iloc[-1]-m[c].iloc[0]:+.4f}, typical sem {(m[c+'_sd']/np.sqrt(nrep)).mean():.4f})")
