"""Fig 2c — top-of-list: per-experiment distribution of each method's top-5%-ranked knockdown.

All violins share ONE density->width scale (not matplotlib's per-violin normalisation), so widths are
directly comparable: every method is scored on the same experiments, each KDE integrates to 1, so a
shared scale means the width at any knockdown level is the true relative density. Each violin is
clipped to its own observed range.
"""
import matplotlib.patheffects as pe
import numpy as np
from scipy.stats import gaussian_kde

from _helpers import ACCENT, BLUE, GREY, INK, TEAL, panel_label

ORDER = ["Random", "OligoWalk·intra", "OligoAI", "TAUSO", "best possible"]
CMAP = {"Random": GREY, "OligoWalk·intra": TEAL, "OligoAI": BLUE, "TAUSO": ACCENT, "best possible": INK}


def draw(ax, D):
    t5 = D.t5
    names = [n for n in ORDER if n in t5]
    cols = [CMAP[n] for n in names]
    xs = np.arange(len(names))
    ygrid = np.linspace(0, 105, 256)
    dens = {n: gaussian_kde(t5[n])(ygrid) for n in names}
    densmax = max(v.max() for v in dens.values())
    for i, (n, c) in enumerate(zip(names, cols)):
        lo, hi = t5[n].min(), t5[n].max()
        m = (ygrid >= lo) & (ygrid <= hi)
        w = dens[n] / densmax * 0.42
        ax.fill_betweenx(ygrid[m], i - w[m], i + w[m], facecolor=c, alpha=0.5, edgecolor=c, linewidth=0.8, zorder=3)
    for i, n in zip(xs, names):
        qa, med, qb = np.percentile(t5[n], [25, 50, 75])
        ax.plot([i, i], [qa, qb], color="#222", lw=5, solid_capstyle="butt", zorder=4)     # IQR
        ax.scatter([i], [med], color="white", s=14, edgecolor="#222", lw=0.5, zorder=5)     # median
        if len(t5[n]) < len(t5["TAUSO"]):
            ax.text(i, 1.5, f"n={len(t5[n])}", ha="center", va="bottom", fontsize=6.5, color="#555", zorder=6)
    # TAUSO median + IQR as light reference lines, values labelled at TAUSO's dot
    halo = [pe.withStroke(linewidth=1.8, foreground="white")]
    ti = names.index("TAUSO")
    tq1, tmed, tq3 = np.percentile(t5["TAUSO"], [25, 50, 75])
    ax.axhline(tmed, color=ACCENT, ls="-", lw=1.0, alpha=0.40, zorder=2)
    for q in (tq1, tq3):
        ax.axhline(q, color=ACCENT, ls="--", lw=0.7, alpha=0.28, zorder=2)
    ax.text(ti + 0.13, tmed, f"{tmed:.1f}", color=INK, fontsize=8.5, va="center", ha="left",
            fontweight="bold", zorder=7, path_effects=halo)
    for q in (tq1, tq3):
        ax.text(ti + 0.13, q, f"{q:.1f}", color="#9aa3b0", fontsize=6.5, va="center", ha="left",
                zorder=7, path_effects=halo)
    ax.set_xticks(xs)
    ax.set_xticklabels([n.replace(" ", "\n").replace("·", "\n·") for n in names], fontsize=8.5)
    ax.set_xlim(-0.6, len(names) - 0.4)
    ax.set_ylabel("top-5% knockdown (%), per experiment", fontsize=10)
    ax.set_ylim(0, 105)
    ax.set_title("Top-of-list: TAUSO's shortlist has a higher floor", fontsize=10.5, fontweight="bold", loc="left", pad=6)
    panel_label(ax, "c")
