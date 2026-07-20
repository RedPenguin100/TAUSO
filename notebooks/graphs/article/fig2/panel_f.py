"""Fig 2f — screening-effort reduction vs random (TAUSO vs OligoAI).

For each screen: to match the knockdown of a method's top-k picks, how many ASOs must a random screen
test (keeping its best k)? fold = N_random / k (OligoAI's own metric). Bars = geometric mean of the
per-screen fold; the dashed line at 1x is no reduction.
"""
import numpy as np

from _helpers import ACCENT, BLUE, GREY, panel_label

ORDER = ["TAUSO", "OligoAI", "random"]
CMAP = {"TAUSO": ACCENT, "OligoAI": BLUE, "random": GREY}


def _geomean(x):
    return float(np.exp(np.mean(np.log(x))))


def draw(ax, D):
    vals = {m: _geomean(D.effort[m].to_numpy()) for m in ORDER}
    xs = np.arange(len(ORDER))
    ax.bar(xs, [vals[m] for m in ORDER], 0.6, color=[CMAP[m] for m in ORDER], edgecolor="white", zorder=3)
    ax.axhline(1.0, color="#888", lw=0.9, ls="--", zorder=2)                 # no-reduction baseline
    for i, m in enumerate(ORDER):
        ax.text(i, vals[m] + 0.06, f"{vals[m]:.1f}×", ha="center", fontsize=11, fontweight="bold", color="#333")
    ax.set_xticks(xs)
    ax.set_xticklabels(ORDER, fontsize=10)
    ax.set_ylabel(f"screening-effort reduction (×)\nto match top-{D.effort_k} picks", fontsize=10)
    ax.set_ylim(0, max(vals.values()) * 1.28)
    ax.set_title("Screening effort — TAUSO saves the most wet-lab work", fontsize=10.5,
                 fontweight="bold", loc="left", pad=6)
    panel_label(ax, "f")
