"""Fig 2d — per-chemistry within-experiment Spearman (TAUSO vs OligoAI)."""
import numpy as np

from _helpers import ACCENT, BLUE, panel_label

CHEMS = ["2'-MOE", "cEt"]


def draw(ax, D):
    x = np.arange(len(CHEMS))
    w = 0.36
    for j, (meth, color) in enumerate([("TAUSO", ACCENT), ("OligoAI", BLUE)]):
        vals = [D.chem_res[c][meth] for c in CHEMS]
        ax.bar(x + (j - 0.5) * w, vals, w, color=color, edgecolor="white", zorder=3, label=meth)
        for i, v in enumerate(vals):
            ax.text(x[i] + (j - 0.5) * w, v + 0.008, f"{v:.2f}", ha="center", fontsize=8.5, color="#333")
    ax.set_xticks(x)
    ax.set_xticklabels(CHEMS, fontsize=11)
    ax.set_ylabel("within-experiment median Spearman", fontsize=10)
    ax.set_title("Per-chemistry ranking — TAUSO's edge is cEt", fontsize=11, fontweight="bold", loc="left", pad=6)
    ax.legend(frameon=False, fontsize=9.5, loc="upper right")
    panel_label(ax, "d")
