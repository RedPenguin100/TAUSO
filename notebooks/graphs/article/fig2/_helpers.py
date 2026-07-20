"""Shared Fig 2 drawing helpers — matplotlib only, no data work."""
import sys
from pathlib import Path

import numpy as np
from matplotlib.ticker import FixedFormatter, FixedLocator

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # notebooks/graphs, for consts
from consts import ACCENT, BLUE, GREY, INK   # noqa: E402  (re-exported for panels)

TEAL = "#2A9D8F"
_XTICKS = [20, 30, 50, 100, 200, 500]


def panel_label(ax, s):
    ax.text(-0.13, 1.05, s, transform=ax.transAxes, fontsize=15, fontweight="bold", va="bottom", color=INK)


def spearman_bar(ax, d, title):
    """Horizontal bar of {method: median Spearman}, TAUSO orange / OligoAI blue / rest grey."""
    items = sorted(d.items(), key=lambda kv: (kv[1] if np.isfinite(kv[1]) else -9))
    labels = [k for k, _ in items]
    vals = [v for _, v in items]
    colors = [ACCENT if k == "TAUSO" else BLUE if k == "OligoAI" else GREY for k in labels]
    yy = np.arange(len(labels))
    ax.barh(yy, vals, color=colors, edgecolor="white", zorder=3)
    ax.set_yticks(yy)
    ax.set_yticklabels(labels, fontsize=9)
    for i, v in enumerate(vals):
        ax.text(v + (0.005 if v >= 0 else -0.005), i, f"{v:.2f}", va="center",
                ha="left" if v >= 0 else "right", fontsize=8, color="#333")
    ax.axvline(0, color="#888", lw=0.8, zorder=2)
    ax.set_xlabel("median Spearman", fontsize=10)
    ax.set_title(title, fontsize=11.5, fontweight="bold", loc="left", pad=6)


def logx_features(ax, kmax):
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator(_XTICKS))
    ax.xaxis.set_major_formatter(FixedFormatter([str(k) for k in _XTICKS]))
    ax.set_xlim(17, kmax * 1.15)
    ax.set_xlabel("number of features (top-K)", fontsize=10)


def droplabels(ax, ks, valof, color):
    for k in ks:
        ax.annotate(f"{valof(k):.2f}", (k, valof(k)), textcoords="offset points", xytext=(0, -12),
                    ha="center", fontsize=7.5, color=color, fontweight="bold")
