"""Helpers for the raw off-target / rRNA feature analysis.

The story: ASOs with off-target burden -- especially binding to ribosomal RNA -- are rarely
effective. Everything is measured WITHIN cohort (custom_id, which nests cell line and is
essentially single-chemistry), because Inhibition(%) is not comparable across experiments.

Features are oriented so HIGHER == MORE off-target burden (the rRNA/RNaseH1 counts are already
non-negative; the off_target_score_* columns are <=0, so they are sign-flipped).
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

INH, COHORT = "Inhibition(%)", "custom_id"
CACHE = Path(__file__).resolve().parent / "data" / "raw_offtarget.parquet"

# label -> (column, sign that makes higher = more burden)
BURDEN = {
    "rRNA (total)":        ("off_target_single_rRNA_total_c1000", +1),
    "rRNA 18S":            ("off_target_single_rRNA_18S_c1000", +1),
    "rRNA 28S":            ("off_target_single_rRNA_28S_c1000", +1),
    "RNase-H1 off-tgt":    ("off_target_single_RNASEH1_c1000", +1),
    "off-target general":  ("off_target_score_general_ARTM_n200_c1000", -1),
    "off-target specific": ("off_target_score_specific_ARTM_n200_c1000", -1),
}
META = [INH, COHORT, "Cell_line", "Modification"]

# plain-language description of each burden feature, for figure captions
BURDEN_DESC = {
    "rRNA (total)":        "binding to ribosomal RNA (18S/28S/5.8S/5S) — the most abundant RNA in the cell",
    "rRNA 18S":            "binding to 18S ribosomal RNA",
    "rRNA 28S":            "binding to 28S ribosomal RNA",
    "RNase-H1 off-tgt":    "binding to the RNase-H1 transcript",
    "off-target general":  "abundance-weighted off-target binding across the transcriptome",
    "off-target specific": "abundance-weighted off-target binding to specifically-expressed transcripts",
}
ACCENT, BLUE, GREY = "#E4572E", "#2E86AB", "#9AA5B1"
SEQ = ["#bdbdbd", "#f3b08a", "#e4572e", "#9e2a0c"]   # none -> high burden (grey -> deep red = worse)


def setup_style():
    plt.rcParams.update({
        "figure.dpi": 120, "font.size": 11, "axes.grid": False,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.titlesize": 13, "axes.titleweight": "bold", "axes.titlelocation": "left",
        "axes.titlepad": 12, "axes.edgecolor": "#bbbbbb",
        "xtick.color": "#555555", "ytick.color": "#555555", "figure.facecolor": "white",
    })


FIGDIR = Path(__file__).resolve().parent / "figures"


def save(fig, name, dpi=200):
    """Save a figure to figures/<name>.png and close it."""
    FIGDIR.mkdir(exist_ok=True)
    path = FIGDIR / f"{name}.png"
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"saved {path}")


def build_cache():
    from notebooks.models.utility import load_and_validate_final_data
    fd, _ = load_and_validate_final_data(version="oligo", load_competition=False)
    cols = META + [c for (c, _s) in BURDEN.values() if c in fd.columns]
    df = fd[cols].copy()
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(CACHE, index=False)
    return df


def load():
    """Slim frame with oriented burden columns, chemistry, and within-cohort residual inhibition."""
    if not CACHE.exists():
        build_cache()
    df = pd.read_parquet(CACHE)
    df["chem"] = df["Modification"].map(lambda m: "cEt" if isinstance(m, str) and "cEt" in m else "2'MOE")
    for name, (col, sgn) in BURDEN.items():
        if col in df.columns:
            df[name] = sgn * df[col].to_numpy("float64")
    df["inh_resid"] = df[INH] - df.groupby(COHORT)[INH].transform("mean")
    return df


def bins(df, name, n=3):
    """Ordered burden categories: 'none' (==0) then n percentile bins of the positive values
    (rank-based, so the heavily-tied rRNA counts still split cleanly)."""
    names = ["low", "high"] if n == 2 else ["low", "mid", "high"] if n == 3 else [f"q{i+1}" for i in range(n)]
    b = df[name].to_numpy("float64")
    out = np.array(["none"] * len(b), dtype=object)
    pos = b > 0
    if pos.sum() >= n:
        g = np.minimum((pd.Series(b[pos]).rank(pct=True).to_numpy() * n).astype(int), n - 1)
        out[pos] = [names[i] for i in g]
    order = ["none"] + [x for x in names if (out == x).any()]
    return pd.Series(pd.Categorical(out, categories=order, ordered=True), index=df.index)


def resid_by_bin(df, name, n=3):
    """Per burden-bin mean within-cohort residual inhibition (+/- SE) and n."""
    g = df.assign(_b=bins(df, name, n)).groupby("_b", observed=True)["inh_resid"]
    return pd.DataFrame({"mean": g.mean(), "se": g.sem(), "n": g.size()})


def burden_log10(df, name):
    """log10(1 + burden): puts the 0-burden ASOs at x=0 and spreads the score's orders of magnitude."""
    return np.log10(1.0 + df[name].to_numpy("float64"))


def binned_resid(df, name, nbins=10, min_n=30):
    """Mean within-cohort residual inhibition vs off-target burden, binned on a log10 axis.

    Equal-width bins of log10(1 + burden); returns one row per non-empty bin with x (bin centre),
    mean, se, n. The smooth alternative to the coarse none/low/mid/high split."""
    b = burden_log10(df, name)
    r = df["inh_resid"].to_numpy("float64")
    edges = np.linspace(b.min(), b.max(), nbins + 1)
    idx = np.clip(np.digitize(b, edges) - 1, 0, nbins - 1)
    rows = []
    for k in range(nbins):
        m = idx == k
        if m.sum() >= min_n:
            rows.append({"x": (edges[k] + edges[k + 1]) / 2, "mean": float(r[m].mean()),
                         "se": float(r[m].std() / np.sqrt(m.sum())), "n": int(m.sum())})
    return pd.DataFrame(rows)


def high_vs_regular(df, name, q=0.25, min_n=12):
    """Per cohort: mean Inhibition(%) of the top-q-burden ASOs ('high') vs the rest ('regular').

    Real inhibition units, but each cohort contributes one high/regular pair, so averaging across
    cohorts controls for the experiment. Returns columns high, reg, n."""
    rows = []
    for _, s in df.groupby(COHORT):
        b = s[name].to_numpy("float64"); y = s[INH].to_numpy("float64")
        if len(y) < min_n or np.ptp(b) == 0:
            continue
        hi = b >= np.quantile(b, 1 - q)
        if hi.sum() < 3 or (~hi).sum() < 3:
            continue
        rows.append((y[hi].mean(), y[~hi].mean(), len(y)))
    return pd.DataFrame(rows, columns=["high", "reg", "n"])


def burden_values(df, name, q=0.25, min_n=12):
    """The burden values that 'high' vs 'regular' actually carry, pooled over the cohorts that
    high_vs_regular uses (so 'high' is the per-cohort top-q). Returns (high, reg) arrays in the
    feature's own units; NaN-safe."""
    hi, reg = [], []
    for _, s in df.groupby(COHORT):
        b = s[name].to_numpy("float64")
        b = b[~np.isnan(b)]
        if len(b) < min_n or np.ptp(b) == 0:
            continue
        m = b >= np.quantile(b, 1 - q)
        if m.sum() < 3 or (~m).sum() < 3:
            continue
        hi.append(b[m]); reg.append(b[~m])
    return np.concatenate(hi), np.concatenate(reg)


def burden_summary(df, name, q=0.25):
    """One-line facts about what counts as 'high' burden for a feature: fraction of all ASOs that
    are nonzero, the median per-cohort 'high' cutoff, and a typical (median) high vs regular value."""
    b = df[name].to_numpy("float64"); b = b[~np.isnan(b)]
    hi, reg = burden_values(df, name, q)
    cuts = [c for c in (
        (lambda x: np.quantile(x, 1 - q) if len(x) >= 12 and np.ptp(x) > 0 else np.nan)
        (s[name].to_numpy("float64")[~np.isnan(s[name].to_numpy("float64"))])
        for _, s in df.groupby(COHORT)) if np.isfinite(c)]
    return {"pct_nonzero": 100 * (b > 0).mean(), "cutoff": float(np.median(cuts)),
            "high": float(np.median(hi)), "reg": float(np.median(reg))}


def cohort_delta(df, name, q=0.25):
    """Per-cohort delta = high_mean - reg_mean (real inhibition %), with the across-cohort summary.

    Returns dict with: delta (per-cohort array), mean & 95% t-CI of the mean delta, Wilcoxon p
    (paired high vs reg), and the means + 95% CIs of high and of regular for dumbbell error bars."""
    from scipy.stats import wilcoxon, t
    d = high_vs_regular(df, name, q)
    delta = (d.high - d.reg).to_numpy("float64")
    n = len(delta)
    ci = lambda v: t.ppf(0.975, len(v) - 1) * v.std(ddof=1) / np.sqrt(len(v))
    return {"delta": delta, "n": n,
            "mean": float(delta.mean()), "ci": float(ci(delta)),
            "p": float(wilcoxon(d.high, d.reg).pvalue),
            "high": float(d.high.mean()), "high_ci": float(ci(d.high.to_numpy())),
            "reg": float(d.reg.mean()),  "reg_ci": float(ci(d.reg.to_numpy()))}


def per_cohort_corr(df, name, min_n=10, min_nonzero=3):
    """Per-cohort Spearman(burden, inhibition); negative => burden flags worse ASOs."""
    out = []
    for _, s in df.groupby(COHORT):
        x = s[name].to_numpy("float64"); y = s[INH].to_numpy("float64")
        if len(y) < min_n or np.ptp(y) == 0 or (x > 0).sum() < min_nonzero or np.ptp(x) == 0:
            continue
        c = spearmanr(x, y)[0]
        if not np.isnan(c):
            out.append(c)
    return np.array(out)


def strong_weak_depletion(df, name, q=0.25):
    """Within each cohort label top-q inhibition 'strong' and bottom-q 'weak'; return the
    fraction of each that carries any of this off-target burden (burden > 0)."""
    d = df.copy()
    d["rank"] = d.groupby(COHORT)[INH].rank(pct=True)
    has = d[name] > 0
    return {"strong": float(has[d["rank"] >= 1 - q].mean()),
            "mid":    float(has[(d["rank"] > q) & (d["rank"] < 1 - q)].mean()),
            "weak":   float(has[d["rank"] <= q].mean())}
