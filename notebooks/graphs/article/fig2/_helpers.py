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


def spearman_bar(ax, d, title, ci=None):
    """Horizontal bar of {method: median Spearman}, TAUSO orange / OligoAI blue / rest grey.

    Method names are left-aligned in the margin (tidy left edge, no clipping); the x-limits carry
    extra room on the negative side so the leftmost bar and its value label clear the axis. `ci`, if
    given, is {method: (lo, hi)} 95% bootstrap intervals drawn as whiskers, with value labels moved
    beyond the interval end.
    """
    items = sorted(d.items(), key=lambda kv: (kv[1] if np.isfinite(kv[1]) else -9))
    labels = [k for k, _ in items]
    vals = np.array([v for _, v in items])
    colors = [ACCENT if k == "TAUSO" else BLUE if k == "OligoAI" else GREY for k in labels]
    yy = np.arange(len(labels))
    ax.barh(yy, vals, color=colors, edgecolor="white", zorder=3)
    ax.set_yticks(yy)
    ax.set_yticklabels(labels, fontsize=9, ha="left")
    ax.tick_params(axis="y", length=0, pad=max(len(l) for l in labels) * 5.4 + 6)
    los = np.array([ci[k][0] for k in labels]) if ci else vals
    his = np.array([ci[k][1] for k in labels]) if ci else vals
    if ci:
        ax.errorbar(vals, yy, xerr=np.vstack([vals - los, his - vals]), fmt="none",
                    ecolor="#444", elinewidth=1.0, capsize=2.5, zorder=4)
    for i, v in enumerate(vals):
        end = his[i] if v >= 0 else los[i]
        ax.text(end + (0.006 if v >= 0 else -0.006), i, f"{v:.2f}", va="center",
                ha="left" if v >= 0 else "right", fontsize=8, color="#333")
    ax.axvline(0, color="#888", lw=0.8, zorder=2)
    lo, hi = min(vals.min(), los.min()), max(vals.max(), his.max())
    span = hi - lo
    ax.set_xlim(lo - (0.16 * span if lo < 0 else 0.02 * span), hi + 0.15 * span)
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
