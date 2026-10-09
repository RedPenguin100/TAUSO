"""What every Figure 2 panel needs: the data in data/, the house style, and the save path.

Holds no panel's content. A panel imports the loader and the style, computes its own numbers and
draws itself, so the panels can be read and rerun one at a time.
"""

from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
OUT = HERE / "out"

INHIBITION, CID, COHORT = "inhibition_percent", "custom_id", "cohort"
MIN_ROWS = 3
"""Smallest screen a within-screen correlation is taken over."""

ACCENT, BLUE, GREY, INK = "#E4572E", "#2E86AB", "#9AA5B1", "#1a1a1a"
TEAL = "#2A9D8F"


def scores():
    """The comparison subset: labels, TAUSO's prediction and every competitor score."""
    path = DATA / "scores.parquet"
    if not path.exists():
        raise SystemExit(f"{path} not found. Run: python {HERE / 'build_data.py'}")
    return pd.read_parquet(path)


def table(name):
    """One of the precomputed sweeps in data/."""
    path = DATA / name
    if not path.exists():
        raise SystemExit(f"{path} not found. See {HERE / 'README.md'}")
    return pd.read_csv(path)


def methods(frame):
    """Every scored method in the frame, TAUSO first."""
    known = [c for c in frame.columns if c not in
             {"index_oligo", INHIBITION, CID, COHORT, "canonical_gene_name", "cell_line", "chemical_pattern"}]
    return ["TAUSO"] + [c for c in known if c != "TAUSO"]


def per_screen_spearman(frame, method, group):
    """Each screen's Spearman between a method's score and measured knockdown."""
    sub = frame[[method, INHIBITION, group]].dropna()
    values = [
        spearmanr(part[method], part[INHIBITION]).correlation
        for _, part in sub.groupby(group)
        if part[method].nunique() > 1 and len(part) >= MIN_ROWS
    ]
    return np.array([v for v in values if np.isfinite(v)])


def median_ci(values, n_boot=2000):
    """(median, lo, hi) with a 95% bootstrap interval, resampling screens rather than ASOs."""
    if not len(values):
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(0)
    draws = np.median(values[rng.integers(0, len(values), size=(n_boot, len(values)))], axis=1)
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return float(np.median(values)), float(lo), float(hi)


def colour_of(method):
    return ACCENT if method == "TAUSO" else BLUE if method == "OligoAI" else GREY


def style():
    """House style, on a transparent canvas so a panel drops into any background."""
    plt.rcParams.update({
        "font.family": "sans-serif", "font.size": 11,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": "#444",
        "figure.facecolor": "none", "axes.facecolor": "none", "savefig.facecolor": "none",
    })


def save(figure, name, caption=None, dpi=300):
    """Write out/<name>.png transparent, a white-background copy beside it, and any caption."""
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"{name}.png"
    figure.savefig(path, dpi=dpi, bbox_inches="tight", transparent=True)
    opaque = Image.new("RGBA", (image := Image.open(path).convert("RGBA")).size, (255, 255, 255, 255))
    Image.alpha_composite(opaque, image).convert("RGB").save(path.with_name(f"{name}_white.png"))
    if caption:
        path.with_suffix(".md").write_text(caption.strip() + "\n")
    plt.close(figure)
    print(f"wrote {path}")


def spearman_bars(ax, medians, intervals, title):
    """Horizontal bars of each method's median Spearman, with bootstrap whiskers."""
    items = sorted(medians.items(), key=lambda kv: kv[1] if np.isfinite(kv[1]) else -9)
    labels = [k for k, _ in items]
    values = np.array([v for _, v in items])
    lows = np.array([intervals[k][0] for k in labels])
    highs = np.array([intervals[k][1] for k in labels])
    y = np.arange(len(labels))

    ax.barh(y, values, color=[colour_of(k) for k in labels], edgecolor="white", zorder=3)
    ax.errorbar(values, y, xerr=np.vstack([values - lows, highs - values]), fmt="none",
                ecolor="#444", elinewidth=1.0, capsize=2.5, zorder=4)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9, ha="left")
    ax.tick_params(axis="y", length=0, pad=max(len(l) for l in labels) * 5.4 + 6)
    for i, value in enumerate(values):
        end = highs[i] if value >= 0 else lows[i]
        ax.text(end + (0.006 if value >= 0 else -0.006), i, f"{value:.2f}", va="center",
                ha="left" if value >= 0 else "right", fontsize=8, color="#333")
    ax.axvline(0, color="#888", lw=0.8, zorder=2)
    low, high = min(values.min(), lows.min()), max(values.max(), highs.max())
    span = high - low
    ax.set_xlim(low - (0.16 * span if low < 0 else 0.02 * span), high + 0.15 * span)
    ax.set_xlabel("median Spearman", fontsize=10)
    ax.set_title(title, fontsize=11.5, fontweight="bold", loc="left", pad=6)
