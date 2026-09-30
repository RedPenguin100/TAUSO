"""TAUSO's parsimony curve against OligoAI on the same held-out TEST split.

Each dot is one rung of the feature descent, scored on the frozen 21,682-row TEST split with
mean +/- sd over 3 seeds. OligoAI is a single horizontal line: its published per-oligo score
put through the same evaluate(), so both sides are measured identically.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ===== TWEAK ME =====
L = Path(__file__).resolve().parent
OUT = L / "tauso_vs_oligoai_parsimony.png"
PANELS = [("test_exp_med", "Within-experiment Spearman\n(median over custom_id)")]
C_TAUSO, C_AI, C_BEST = "#1f5f8b", "#7a4a21", "#0d3350"
C_ORANGE, C_PURPLE = "#e8830c", "#7b4f9e"
MODEL_STYLE = {"MED clean_exp": (C_ORANGE, "--"), "MED reg": (C_PURPLE, "-")}
SMOOTH_SIGMA = 0.04  # gaussian width in ln(n), i.e. about 4% of the feature count
FIGSIZE = (7.5, 4.8)
# =====================

te = pd.read_csv(L / "test_curve.csv")
other = pd.read_csv(L / "other_models_test.csv")
ai = json.loads((L / "oligoai_test_panel.json").read_text())

fig, axes = plt.subplots(1, len(PANELS), figsize=FIGSIZE, squeeze=False)
axes = axes[0]
for ax, (col, label) in zip(axes, PANELS):
    key = col.replace("test_", "")
    ln = np.log(te.n.to_numpy(float))
    w = np.exp(-0.5 * ((ln[:, None] - ln[None, :]) / SMOOTH_SIGMA) ** 2)
    smooth = (w @ te[col].to_numpy()) / w.sum(1)
    ax.plot(te.n, te[col], "o", ms=2.6, color=C_TAUSO, alpha=0.3, lw=0, label="TAUSO, each feature count")
    ax.plot(te.n, smooth, "-", lw=2.2, color=C_TAUSO, label="TAUSO, smoothed")
    ax.axhline(ai[key], color=C_AI, ls="-", lw=1.6, label=f"OligoAI = {ai[key]:.4f}")
    if key == "exp_med":
        full = te.loc[te.n == te.n.max()].iloc[0]
        ax.axhline(full[col], color=C_ORANGE, ls="-", lw=1.6,
                   label=f"all {int(full.n)} features = {full[col]:.4f} (3 seeds)")
        for o in other[other.label.isin(MODEL_STYLE)].itertuples():
            color, ls = MODEL_STYLE[o.label]
            ax.axhline(o.test_exp_med, color=color, ls=ls, lw=1.6,
                       label=f"{o.label} = {o.test_exp_med:.4f} ({o.seeds} seeds)")
    best = te.loc[te[col].idxmax()]
    ax.plot(best.n, best[col], marker="*", ms=15, color=C_BEST, lw=0, zorder=5,
            label=f"best so far: n={int(best.n)}, {best[col]:.4f}")
    ax.set_xscale("log")
    ax.invert_xaxis()
    ticks = [t for t in (679, 500, 400, 300, 200, 150, 100, 60, 40, 25, 15, 10)
             if te.n.min() * 0.95 <= t <= te.n.max() * 1.05]
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(t) for t in ticks])
    ax.minorticks_off()
    ax.set_xlabel("features in the model (descending)")
    ax.set_ylabel(label, fontsize=9)
    lo = min(te[col].min(), ai[key])
    hi = max(te[col].max(), ai[key])
    pad = 0.12 * (hi - lo)
    ax.set_ylim(lo - pad, hi + 1.6 * pad)
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(fontsize=8, loc="center left", framealpha=0.9)

fig.suptitle(f"TAUSO feature descent vs OligoAI, held-out TEST (n={ai['n_rows']:,} oligos)",
             fontsize=11)
fig.tight_layout()
fig.savefig(OUT, dpi=200)
print(f"wrote {OUT}")
for ax, (col, _) in zip(axes, PANELS):
    key = col.replace("test_", "")
    print(f"{col}: TAUSO {te[col].min():.4f}-{te[col].max():.4f} over n={te.n.max()}-{te.n.min()}"
          f"   OligoAI {ai[key]:.4f}   margin at best {te[col].max()-ai[key]:+.4f}")
