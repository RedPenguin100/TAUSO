"""Fig 2f — screening-effort reduction vs random (TAUSO vs OligoAI).

For each screen: to match the knockdown of a method's top-k picks, how many ASOs must a random screen
test (keeping its best k)? fold = N_random / k (OligoAI's own metric). Bars = geometric mean of the
per-screen fold; whiskers = 95% bootstrap CI resampling screens. The dashed line at 1x is no reduction;
random sits just above it (the fold is floored at 1x, so its geomean is slightly inflated).
"""
import numpy as np

from _helpers import ACCENT, BLUE, GREY, panel_label

ORDER = ["TAUSO", "OligoAI", "random"]
CMAP = {"TAUSO": ACCENT, "OligoAI": BLUE, "random": GREY}


def draw(ax, D):
    stat = D.effort_stat                                    # method -> (geomean, ci_lo, ci_hi)
    xs = np.arange(len(ORDER))
    gm = [stat[m][0] for m in ORDER]
    ax.bar(xs, gm, 0.6, color=[CMAP[m] for m in ORDER], edgecolor="white", zorder=3)
    err = np.array([[stat[m][0] - stat[m][1], stat[m][2] - stat[m][0]] for m in ORDER]).T
    ax.errorbar(xs, gm, yerr=err, fmt="none", ecolor="#444", elinewidth=1.1, capsize=4, zorder=4)
    ax.axhline(1.0, color="#888", lw=0.9, ls="--", zorder=2)                 # no-reduction baseline
    top = max(stat[m][2] for m in ORDER)
    for i, m in enumerate(ORDER):
        ax.text(i, stat[m][2] + 0.03 * top, f"{stat[m][0]:.1f}×", ha="center",
                fontsize=11, fontweight="bold", color="#333")
    ax.set_xticks(xs)
    ax.set_xticklabels(ORDER, fontsize=10)
    ax.set_ylabel(f"screening-effort reduction (×)\nto match top-{D.effort_k} picks", fontsize=10)
    ax.set_ylim(0, top * 1.16)
    ax.set_title("Screening effort — TAUSO saves the most wet-lab work", fontsize=10.5,
                 fontweight="bold", loc="left", pad=6)
    panel_label(ax, "f")
