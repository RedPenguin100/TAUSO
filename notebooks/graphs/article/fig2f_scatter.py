"""Fig 2 panel f (standalone preview): paired per-screen ranking, TAUSO vs OligoAI.

One dot per screen (custom_id) on the Fig 2 comparison subset: x = OligoAI's within-screen Spearman
vs measured knockdown, y = TAUSO's. Points above the diagonal are screens where TAUSO ranks better.

    python fig2f_scatter.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, wilcoxon

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "models"))
from consts import save, style, ACCENT, BLUE, INK   # noqa: E402
import _modeldata as M                              # noqa: E402

CAT = "article"


def per_screen(d, sign):
    rows = []
    for cid, s in d.groupby(M.CID):
        def sp(col, sg=1.0):                                        # same rule as _modeldata.med (a/b/d/e)
            x = s[[col, M.INH]].dropna()
            return spearmanr(x[col], x[M.INH]).correlation * sg if x[col].nunique() > 1 and len(x) >= 3 else np.nan
        t = sp("tauso")
        o = sp("oligo_ai_score", sign["OligoAI"])
        if np.isfinite(t) and np.isfinite(o):
            rows.append({"tauso": t, "oligoai": o, "n": len(s)})
    return pd.DataFrame(rows)


def main():
    d, present = M.load()
    sign = M.orient(d, present)
    p = per_screen(d, sign)
    win = (p.tauso > p.oligoai).mean()
    delta = p.tauso - p.oligoai
    stat, pval = wilcoxon(delta)

    style()
    fig, ax = plt.subplots(figsize=(5.4, 5.4))
    lim = (-0.65, 1.0)
    ax.plot(lim, lim, color="#888", lw=1.0, ls="--", zorder=1)                 # y = x
    above = p.tauso > p.oligoai
    ax.scatter(p.oligoai[above], p.tauso[above], s=16, color=ACCENT, alpha=0.6,
               edgecolor="none", zorder=3, label=f"TAUSO better ({100*win:.0f}%)")
    ax.scatter(p.oligoai[~above], p.tauso[~above], s=16, color=BLUE, alpha=0.55,
               edgecolor="none", zorder=3, label=f"OligoAI better ({100*(1-win):.0f}%)")
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal")
    ax.set_xlabel("OligoAI — within-screen Spearman", fontsize=10)
    ax.set_ylabel("TAUSO — within-screen Spearman", fontsize=10)
    ax.set_title("Per-screen ranking: TAUSO vs OligoAI", fontsize=11.5, fontweight="bold", loc="left", pad=6)
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    ax.text(0.03, 0.97, f"n = {len(p)} screens\nmedian Δ = {np.median(delta):+.2f}\n"
            f"Wilcoxon p = {pval:.1e}", transform=ax.transAxes, va="top", ha="left",
            fontsize=8.5, color=INK)

    cap = (
        f"Per-screen ranking, TAUSO vs OligoAI (n = {len(p)} screens on the Fig 2 comparison subset). "
        f"Each dot is one screen (custom_id); axes are the within-screen Spearman of each method's score "
        f"vs measured knockdown. Points above the diagonal are screens TAUSO ranks better: {100*win:.0f}% "
        f"of screens (median advantage {np.median(delta):+.2f}, paired Wilcoxon p = {pval:.1e})."
    )
    save(fig, CAT, "fig2f_scatter", caption=cap)
    print(f"fig2f built | n={len(p)} screens | TAUSO wins {100*win:.0f}% | median dcid {np.median(delta):+.3f} | p={pval:.2e}")


if __name__ == "__main__":
    main()
