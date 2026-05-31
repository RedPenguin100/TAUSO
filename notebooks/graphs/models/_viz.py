"""Shared style + plotting + metric helpers for the model-comparison notebooks.

Thin, reusable pieces so the notebooks stay short. All compute functions take the test
predictions table written by Evaluation/evaluate.py (one row per test oligo, one column per
scorer) and aggregate per patent cohort (custom_id).
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, rankdata

INHIB, COHORT, MODEL, OLIGOAI = "Inhibition(%)", "custom_id", "TAUSO", "oligo_ai_score"
RESULTS = Path(__file__).resolve().parents[2] / "models/Evaluation/results"

ACCENT, BLUE, GREY = "#E4572E", "#2E86AB", "#9AA5B1"
DISPLAY = {
    MODEL: "TAUSO", OLIGOAI: "OligoAI", "PFRED_PLS": "PFRED", "PFRED_SVM": "PFRED·SVM",
    "OW_Overall": "OligoWalk", "OW_Tm": "OligoWalk·Tm", "OW_Intra_Oligo": "OligoWalk·intra",
    "OW_Duplex": "OligoWalk·duplex", "sfold_accessibility": "sfold",
    "miranda_score": "miRanda", "miranda_energy": "miRanda·E", "cohort_mean": "cohort-mean",
}
HEADLINE = [MODEL, OLIGOAI, "OW_Intra_Oligo", "sfold_accessibility", "PFRED_PLS", "miranda_energy"]

# multi-output tools: pick each tool's best-ranking column for the headline figures
FAMILIES = {"PFRED": ["PFRED_PLS", "PFRED_SVM"],
            "OligoWalk": ["OW_Overall", "OW_Tm", "OW_Intra_Oligo", "OW_Duplex"],
            "miRanda": ["miranda_score", "miranda_energy"]}
SINGLES = [MODEL, OLIGOAI, "sfold_accessibility"]


def disp(s):  return DISPLAY.get(s, s)
def color(s): return ACCENT if s == MODEL else (BLUE if s == OLIGOAI else GREY)


def setup_style():
    plt.rcParams.update({
        "figure.dpi": 120, "font.size": 11, "axes.grid": False,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.titlesize": 13, "axes.titleweight": "bold", "axes.titlelocation": "left",
        "axes.titlepad": 12, "axes.edgecolor": "#bbbbbb",
        "xtick.color": "#555555", "ytick.color": "#555555", "figure.facecolor": "white",
    })


def load_results(d=RESULTS):
    d = Path(d)
    return {
        "summary":   pd.read_csv(d / "metrics_summary.csv"),
        "per_group": pd.read_csv(d / "per_group_spearman.csv"),
        "breakdown": pd.read_csv(d / "breakdown_metrics.csv"),
        "preds":     pd.read_parquet(d / "test_predictions.parquet"),
    }


def order(df, metric, subset=None):
    """Scorers present for `metric`, TAUSO first then descending."""
    s = df.dropna(subset=[metric]).set_index("scorer")[metric]
    if subset is not None:
        s = s[[x for x in subset if x in s.index]]
    rest = [x for x in s.sort_values(ascending=False).index if x != MODEL]
    return ([MODEL] + rest) if MODEL in s.index else rest


def hbar(ax, scorers, values, title="", xlabel="", fmt="{:.2f}"):
    """Horizontal bar chart, TAUSO accented, value labels, bold left title."""
    y = np.arange(len(scorers))[::-1]
    ax.barh(y, values, color=[color(s) for s in scorers])
    ax.set_yticks(y); ax.set_yticklabels([disp(s) for s in scorers])
    for yi, v in zip(y, values):
        if np.isfinite(v):
            ax.text(v, yi, ("  " + fmt.format(v)) if v >= 0 else (fmt.format(v) + "  "),
                    va="center", ha="left" if v >= 0 else "right", fontsize=9, color="#333")
    ax.set_title(title); ax.set_xlabel(xlabel); ax.axvline(0, color="#cccccc", lw=0.8)
    lo, hi = float(np.nanmin(values)), float(np.nanmax(values))
    rng = hi - min(0, lo)
    ax.set_xlim(min(0, lo) - 0.16 * rng - 0.01, hi + 0.16 * rng + 0.01)
    ax.margins(y=0.02)
    return ax


# ----------------------------------------------------------------- per-cohort compute

def _cohort_arrays(preds, scorer, cohort=COHORT, min_n=8):
    for _, sub in preds.groupby(cohort):
        y = sub[INHIB].to_numpy("float64"); s = sub[scorer].to_numpy("float64")
        ok = np.isfinite(y) & np.isfinite(s)
        if ok.sum() >= min_n:
            yield y[ok], s[ok]


def add_cohort_mean(preds, cohort=COHORT):
    """Add a 'cohort_mean' column: each oligo scored by its cohort's mean inhibition
    (the leakage baseline — high global Spearman, zero ranking signal within a cohort)."""
    out = preds.copy()
    out["cohort_mean"] = out.groupby(cohort)[INHIB].transform("mean")
    return out


def global_spearman(preds, scorer):
    y = preds[INHIB].to_numpy("float64"); p = preds[scorer].to_numpy("float64")
    ok = np.isfinite(y) & np.isfinite(p)
    return float(spearmanr(y[ok], p[ok])[0])


def hit_rate_curve(preds, scorer, threshold=70, ks=(1, 3, 5, 10, 15, 20), cohort=COHORT):
    """Mean (over cohorts) precision and enrichment of the top-k by `scorer`; a hit is
    inhibition >= threshold. Enrichment = top-k hit-rate / cohort base hit-rate."""
    rows = []
    for y, s in _cohort_arrays(preds, scorer, cohort):
        base = (y >= threshold).mean()
        if base == 0:
            continue
        o = np.argsort(-s)
        for k in ks:
            kk = min(k, len(y))
            prec = (y[o[:kk]] >= threshold).mean()
            rows.append((k, prec, prec / base))
    g = pd.DataFrame(rows, columns=["k", "precision", "enrichment"])
    return g.groupby("k").mean().reset_index()


def base_hit_rate(preds, threshold=70, cohort=COHORT):
    rs = [(y >= threshold).mean() for y, _ in _cohort_arrays(preds, INHIB, cohort)]
    return float(np.mean(rs)) if rs else float("nan")


def gap_closed(preds, scorer, k=10, cohort=COHORT):
    """Per cohort: (top-k potency - random) / (oracle top-k - random). 1 = oracle, 0 = random."""
    out = []
    for y, s in _cohort_arrays(preds, scorer, cohort, min_n=k + 1):
        rand, oracle = y.mean(), np.sort(y)[-k:].mean()
        if oracle - rand < 1.0:
            continue
        out.append((y[np.argsort(-s)[:k]].mean() - rand) / (oracle - rand))
    return np.array(out)


def best_rank_cdf(preds, scorer, cohort=COHORT):
    """Predicted rank (1 = top) of each cohort's truly most-potent oligo.

    Tie-aware (average rank): a scorer that gives many oligos the same value gets no credit for
    'ranking' the best one first -- it gets the expected rank over the tied group. Without this,
    heavily-tied biophysical scores (e.g. OligoWalk intra-oligo, ~50 oligos tied per cohort) look
    spuriously perfect."""
    ranks = [rankdata(-s, method="average")[np.argmax(y)]
             for y, s in _cohort_arrays(preds, scorer, cohort)]
    return np.array(ranks)


def per_cohort_spearman(preds, scorer, cohort=COHORT, min_n=8):
    """Per-cohort Spearman values of one scorer (for distribution overlays)."""
    out = []
    for y, s in _cohort_arrays(preds, scorer, cohort, min_n):
        if np.ptp(s) > 0:
            c = spearmanr(y, s)[0]
            if np.isfinite(c):
                out.append(c)
    return np.array(out)


def best_representatives(summary, grouping="patent", metric="spearman_median"):
    """Each multi-output tool's best-ranking column + a note (tool -> (column, n_variants))."""
    d = summary[summary.grouping == grouping].set_index("scorer")[metric]
    reps, note = [], {}
    for tool, cols in FAMILIES.items():
        present = [c for c in cols if c in d.index and np.isfinite(d[c])]
        if present:
            best = max(present, key=lambda c: d[c])
            reps.append(best)
            note[tool] = (best, len(present))
    return reps, note


def rep_caption(note):
    return ("multi-output tools shown by their best-ranking column: "
            + ", ".join(f"{disp(b)} (best of {n})" for _, (b, n) in note.items()) + ".")


def hbar_dist(ax, scorers, medians, dists, title="", xlabel="", seed=0, colors=None):
    """Median bar + per-cohort distribution (jittered dots), labelled with the median.

    Rows are labelled by `scorers`; bar colour is `color(scorer)` unless `colors` (one per row)
    is given (e.g. a single accent when every row is the same model across cell lines)."""
    rng = np.random.default_rng(seed)
    cols = colors if colors is not None else [color(s) for s in scorers]
    y = np.arange(len(scorers))[::-1]
    for yi, c, med, d in zip(y, cols, medians, dists):
        ax.barh(yi, med, height=0.6, color=c, alpha=0.45, zorder=1)
        if len(d):
            ax.scatter(d, yi + (rng.random(len(d)) - 0.5) * 0.5, s=7,
                       color=c, alpha=0.30, edgecolors="none", zorder=2)
        ax.text(med, yi, f"  {med:.2f}", fontsize=8.5, va="center", ha="left",
                color="#222", zorder=6,
                bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.55))
    ax.set_yticks(y); ax.set_yticklabels([disp(s) for s in scorers])
    ax.set_title(title); ax.set_xlabel(xlabel); ax.axvline(0, color="#cccccc", lw=0.8)
    allv = np.concatenate([d for d in dists if len(d)] + [np.asarray(medians, float)])
    lo, hi = float(np.nanmin(allv)), float(np.nanmax(allv))
    rng2 = hi - min(0, lo)
    ax.set_xlim(min(0, lo) - 0.10 * rng2 - 0.02, hi + 0.12 * rng2 + 0.02)
    return ax


def reliability(preds, scorer, bins=10, min_bin=20):
    """Calibration curve: mean predicted vs mean observed inhibition per predicted-decile."""
    y = preds[INHIB].to_numpy("float64"); p = preds[scorer].to_numpy("float64")
    ok = np.isfinite(y) & np.isfinite(p); y, p = y[ok], p[ok]
    idx = np.clip(np.digitize(p, np.linspace(0, 100, bins + 1)) - 1, 0, bins - 1)
    rows = [(p[idx == b].mean(), y[idx == b].mean(), int((idx == b).sum()))
            for b in range(bins) if (idx == b).sum() >= min_bin]
    return pd.DataFrame(rows, columns=["pred", "obs", "n"])
