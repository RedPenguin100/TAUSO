"""Fig 2 — TAUSO vs the field (performance mega-panel, DRAFT).

a) within-experiment (custom_id) Spearman     b) within-cohort (gene×cell) Spearman
c) top-of-list: median over experiments of each method's top-5%-ranked knockdown, random floor -> oracle ceiling
d) per-chemistry (2'-MOE / cEt)
e) parsimony — cross-validation (selection)   f) parsimony — held-out test (confirmation vs competitors)
Held-out test, strict 5-10-5 MOE / 3-10-3 cEt gapmers, OligoAI-common rows. TAUSO = held-out model.
Run: python fig2_performance.py
"""
import os
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter
from scipy.stats import gaussian_kde, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))                 # graphs/ for consts
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "models"))      # _modeldata
from consts import save, style, ACCENT, BLUE, GREY, INK
import _modeldata as M
from notebooks.models.evaluate import evaluate
from notebooks.data.OligoAI.parse_chemistry import assign_chemistry
from tauso.data.consts import CHEMICAL_PATTERN

CAT = "article"
WITH_PARSIMONY = False        # a-d only when False; a-f (adds parsimony e,f) when True
AB_IQR = bool(os.environ.get("FIG2_AB_IQR"))   # add IQR whiskers to ranking panels a/b (saved as *_abIQR)
INK_ = INK
INH, CID = M.INH, M.CID
TEAL = "#2A9D8F"
GRAPHS = Path(__file__).resolve().parents[1]
PARS_TEST = GRAPHS / "models/parsimony_test_metrics.csv"
CV_TXT = GRAPHS.parents[0] / "models/results/feature_selection_clean_exp.txt"
LOWK, XTICKS = [20, 30, 40, 50], [20, 30, 50, 100, 200, 500]


def chem_of(p):
    p = str(p)
    if "M" in p and "C" in p: return "mixmer"
    if "M" in p: return "2'-MOE"
    if "C" in p: return "cEt"
    if len(p) > 0 and set(p) <= set("d"): return "PS-DNA"
    return "other"


# ---- load all competitors + TAUSO + chemistry on the held-out test split ----
avg = pd.read_csv(M.AVG, usecols=[M.IDX, INH, CID, "aso_sequence", "split", CHEMICAL_PATTERN, M.GENE, M.CELL])
avg["chem"] = avg[CHEMICAL_PATTERN].map(chem_of)
avg["is_strict"] = avg[CHEMICAL_PATTERN].isin(M.STRICT)
df = avg
present = {}
for lbl, stem in M.COMP.items():
    try:
        df = df.merge(M._shard(stem), on=M.IDX, how="left"); present[lbl] = stem
    except FileNotFoundError:
        pass
df = df.merge(pd.read_parquet(M.TAUSO_PRED), on=M.IDX, how="left"); present["TAUSO"] = "tauso"
df = df[(df["split"] == "test") & df[INH].notna()].reset_index(drop=True)
df["gxc"] = df[M.GENE].astype(str) + "|" + df[M.CELL].astype(str)
S = df[df["is_strict"] & df["oligo_ai_score"].notna()].reset_index(drop=True)
sign = M.orient(S, present)
cid = {lbl: M.med(S, col, CID) * sign[lbl] for lbl, col in present.items()}
gxc = {lbl: M.med(S, col, "gxc") * sign[lbl] for lbl, col in present.items()}
n_cid, n_gxc = S[CID].nunique(), S["gxc"].nunique()


def _grp_iqr(score, grp, s):
    """(Q1, median, Q3) of the per-group (custom_id / gene×cell) Spearman of `score` vs inhibition, sign-oriented."""
    d = S[[score, INH, grp]].dropna()
    rs = [spearmanr(g[score], g[INH]).correlation for _, g in d.groupby(grp)
          if g[score].nunique() > 1 and len(g) >= 3]
    rs = np.array([r for r in rs if np.isfinite(r)]) * s
    return tuple(np.percentile(rs, [25, 50, 75])) if len(rs) else (np.nan, np.nan, np.nan)


cid_iqr = {lbl: _grp_iqr(col, CID, sign[lbl]) for lbl, col in present.items()} if AB_IQR else None
gxc_iqr = {lbl: _grp_iqr(col, "gxc", sign[lbl]) for lbl, col in present.items()} if AB_IQR else None


def top5_dist(score, s=1.0):
    out = []
    for _, g in S.groupby(CID):
        n = len(g)
        if n >= 10 and g[INH].nunique() > 1:
            k = max(1, int(round(0.05 * n)))
            if score == "_random": out.append(g[INH].mean())
            elif score == "_oracle": out.append(g.nlargest(k, INH)[INH].mean())
            else:                                              # top-k by sign-oriented score; skip experiments with no score
                sel = (g[score] * s).dropna().nlargest(k)
                if len(sel): out.append(g.loc[sel.index, INH].mean())
    return np.array(out)


_t5spec = [("Random", "_random", 1.0),
           ("OligoWalk·intra", present.get("OligoWalk·intra"), sign.get("OligoWalk·intra", 1.0)),
           ("OligoAI", "oligo_ai_score", sign["OligoAI"]), ("TAUSO", "tauso", sign["TAUSO"]),
           ("best possible", "_oracle", 1.0)]
t5 = {n: top5_dist(sc, sg) for n, sc, sg in _t5spec if sc is not None}
chem_res = {c: {"TAUSO": M.med(S[S.chem == c], "tauso", CID), "OligoAI": M.med(S[S.chem == c], "oligo_ai_score", CID)}
            for c in ["2'-MOE", "cEt"]}

# ---- parsimony data ----
if WITH_PARSIMONY:
    cv = {}
    for line in CV_TXT.read_text().splitlines():
        mt = re.match(r"K=(\d+)\s+([\d.]+)±([\d.]+)", line)
        if mt:
            cv[int(mt.group(1))] = (float(mt.group(2)), float(mt.group(3)))
    cvK = np.array(sorted(cv)); cvm = np.array([cv[k][0] for k in cvK]); cvs = np.array([cv[k][1] for k in cvK])
    pt = pd.read_csv(PARS_TEST).sort_values("K"); ptmap = dict(zip(pt["K"], pt["exp_med"]))
    ptK, ptm = pt["K"].to_numpy(), pt["exp_med"].to_numpy()
    _y, _c, _g = S[INH].to_numpy(np.float64), S[CID].to_numpy(), S["gxc"].to_numpy()
    _em = lambda stem, s: evaluate((s * S[stem]).to_numpy(np.float64), _y, _c, _g)["exp_med"]
    oai_em = _em("oligo_ai_score", sign["OligoAI"])
    _ow = {lbl: _em(present[lbl], sign[lbl]) for lbl in present if "OligoWalk" in lbl or lbl == "sfold"}
    ow_lbl = max(_ow, key=_ow.get); ow_em = _ow[ow_lbl]


def lbl(ax, s):
    ax.text(-0.13, 1.05, s, transform=ax.transAxes, fontsize=15, fontweight="bold", va="bottom", color=INK)


def spearman_panel(ax, d, title, iqr=None):
    items = sorted(d.items(), key=lambda kv: (kv[1] if np.isfinite(kv[1]) else -9))
    labels = [k for k, _ in items]; vals = [v for _, v in items]
    colors = [ACCENT if k == "TAUSO" else BLUE if k == "OligoAI" else GREY for k in labels]
    yy = np.arange(len(labels))
    ax.barh(yy, vals, color=colors, edgecolor="white", zorder=3)
    if iqr is not None:                                          # IQR whisker (Q1-Q3) + median tick per method
        for i, k in enumerate(labels):
            q1, md, q3 = iqr[k]
            ax.plot([q1, q3], [i, i], color="#2a2a2a", lw=1.2, alpha=0.65, solid_capstyle="butt", zorder=4)
            for q in (q1, q3):
                ax.plot([q, q], [i - 0.18, i + 0.18], color="#2a2a2a", lw=1.2, alpha=0.65, zorder=4)
            ax.scatter([md], [i], s=9, color="white", edgecolor="#2a2a2a", lw=0.5, zorder=5)
    ax.set_yticks(yy); ax.set_yticklabels(labels, fontsize=9)
    for i, v in enumerate(vals):
        if iqr is not None:
            ax.text(iqr[labels[i]][2] + 0.012, i, f"{v:.2f}", va="center", ha="left", fontsize=7.5, color="#333", zorder=6)
        else:
            ax.text(v + (0.005 if v >= 0 else -0.005), i, f"{v:.2f}", va="center",
                    ha="left" if v >= 0 else "right", fontsize=8, color="#333")
    ax.axvline(0, color="#888", lw=0.8, zorder=2)
    if iqr is not None:                                          # widen x so whiskers + labels fit
        _q1s = [iqr[k][0] for k in labels]; _q3s = [iqr[k][2] for k in labels]
        ax.set_xlim(min(_q1s) - 0.04, max(_q3s) + 0.12)
    ax.set_xlabel("median Spearman", fontsize=10)
    ax.set_title(title, fontsize=11.5, fontweight="bold", loc="left", pad=6)


def logx(ax):
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator(XTICKS)); ax.xaxis.set_major_formatter(FixedFormatter([str(k) for k in XTICKS]))
    ax.set_xlim(17, 560); ax.set_xlabel("number of features (top-K)", fontsize=10)


def droplabels(ax, valof, color):
    for k in LOWK:
        ax.annotate(f"{valof(k):.2f}", (k, valof(k)), textcoords="offset points", xytext=(0, -12),
                    ha="center", fontsize=7.5, color=color, fontweight="bold")


style()
fig = plt.figure(figsize=(13.4, 14.6) if WITH_PARSIMONY else (13.4, 9.6))
gs = fig.add_gridspec(3 if WITH_PARSIMONY else 2, 2)

axa = fig.add_subplot(gs[0, 0]); spearman_panel(axa, cid, "Within-experiment ranking (per patent / custom_id)", iqr=cid_iqr); lbl(axa, "a")
axb = fig.add_subplot(gs[0, 1]); spearman_panel(axb, gxc, "Within-cohort ranking (per gene × cell-line)", iqr=gxc_iqr); lbl(axb, "b")

# c — top-of-list: full per-experiment distribution (violin) of each method's top-5%-ranked knockdown.
# All four violins share ONE density->width scale (not matplotlib's per-violin max-width normalisation), so
# widths are directly comparable = honest proportions: every method is scored on the SAME experiments and each
# KDE integrates to 1, so a shared scale means equal area and the width at any knockdown level is the true
# relative density there. Each violin is clipped to its own observed range.
axc = fig.add_subplot(gs[1, 0])
_order = ["Random", "OligoWalk·intra", "OligoAI", "TAUSO", "best possible"]
_cmap = {"Random": GREY, "OligoWalk·intra": TEAL, "OligoAI": BLUE, "TAUSO": ACCENT, "best possible": INK}
names = [n for n in _order if n in t5]; cols = [_cmap[n] for n in names]
xs = np.arange(len(names))
ygrid = np.linspace(0, 105, 256)
dens = {n: gaussian_kde(t5[n])(ygrid) for n in names}
kmax = max(d.max() for d in dens.values())
for i, (n, c) in enumerate(zip(names, cols)):
    lo, hi = t5[n].min(), t5[n].max()
    m = (ygrid >= lo) & (ygrid <= hi)
    w = dens[n] / kmax * 0.42
    axc.fill_betweenx(ygrid[m], i - w[m], i + w[m], facecolor=c, alpha=0.5, edgecolor=c, linewidth=0.8, zorder=3)
for i, n in zip(xs, names):
    qa, med, qb = np.percentile(t5[n], [25, 50, 75])
    axc.plot([i, i], [qa, qb], color="#222", lw=5, solid_capstyle="butt", zorder=4)      # IQR
    axc.scatter([i], [med], color="white", s=14, edgecolor="#222", lw=0.5, zorder=5)      # median
    if len(t5[n]) < len(t5["TAUSO"]):                                                     # note reduced coverage
        axc.text(i, 1.5, f"n={len(t5[n])}", ha="center", va="bottom", fontsize=6.5, color="#555", zorder=6)
# TAUSO median + IQR as light horizontal reference lines across the panel (so competitors visibly sit below),
# with the values labelled at TAUSO's own dot: median bold, IQR edges faint.
import matplotlib.patheffects as _pe
_halo = [_pe.withStroke(linewidth=1.8, foreground="white")]
_ti = names.index("TAUSO")
_tq1, _tmed, _tq3 = np.percentile(t5["TAUSO"], [25, 50, 75])
axc.axhline(_tmed, color=ACCENT, ls="-", lw=1.0, alpha=0.40, zorder=2)                      # TAUSO median
for _q in (_tq1, _tq3):
    axc.axhline(_q, color=ACCENT, ls="--", lw=0.7, alpha=0.28, zorder=2)                    # TAUSO IQR edges
axc.text(_ti + 0.13, _tmed, f"{_tmed:.1f}", color=INK_, fontsize=8.5, va="center", ha="left",
         fontweight="bold", zorder=7, path_effects=_halo)                                   # median value at the dot
for _q in (_tq1, _tq3):
    axc.text(_ti + 0.13, _q, f"{_q:.1f}", color="#9aa3b0", fontsize=6.5, va="center", ha="left",
             zorder=7, path_effects=_halo)                                                  # faint IQR values
axc.set_xticks(xs); axc.set_xticklabels([n.replace(" ", "\n").replace("·", "\n·") for n in names], fontsize=8.5)
axc.set_xlim(-0.6, len(names) - 0.4)
axc.set_ylabel("top-5% knockdown (%), per experiment", fontsize=10); axc.set_ylim(0, 105)
axc.set_title("Top-of-list: TAUSO's shortlist has a higher floor", fontsize=10.5, fontweight="bold", loc="left", pad=6)
lbl(axc, "c")

# d — per-chemistry
axd = fig.add_subplot(gs[1, 1])
chems = ["2'-MOE", "cEt"]; x = np.arange(2); w = 0.36
for j, (meth, color) in enumerate([("TAUSO", ACCENT), ("OligoAI", BLUE)]):
    vals = [chem_res[c][meth] for c in chems]
    axd.bar(x + (j - 0.5) * w, vals, w, color=color, edgecolor="white", zorder=3, label=meth)
    for i, v in enumerate(vals):
        axd.text(x[i] + (j - 0.5) * w, v + 0.008, f"{v:.2f}", ha="center", fontsize=8.5, color="#333")
axd.set_xticks(x); axd.set_xticklabels(chems, fontsize=11)
axd.set_ylabel("within-experiment median Spearman", fontsize=10)
axd.set_title("Per-chemistry ranking — TAUSO's edge is cEt", fontsize=11, fontweight="bold", loc="left", pad=6)
axd.legend(frameon=False, fontsize=9.5, loc="upper right")
lbl(axd, "d")

if WITH_PARSIMONY:
    # e — parsimony, cross-validation (selection)
    axe = fig.add_subplot(gs[2, 0])
    axe.axvline(450, color="#ccc", lw=1.0, ls=":", zorder=1)
    axe.fill_between(cvK, cvm - cvs, cvm + cvs, color="#5B6B7B", alpha=0.16, lw=0, zorder=2)
    axe.plot(cvK, cvm, color="#5B6B7B", lw=2.0, zorder=3)
    axe.scatter([450], [cv[450][0]], color=ACCENT, s=170, marker="*", edgecolor=INK, lw=0.8, zorder=5)
    axe.annotate("chosen K=450", (450, cv[450][0]), textcoords="offset points", xytext=(-6, 9), ha="right", fontsize=8.5, color=ACCENT, fontweight="bold")
    droplabels(axe, lambda k: cv[k][0], "#5B6B7B")
    logx(axe); axe.set_ylabel("within-experiment Spearman", fontsize=10)
    axe.set_title("Model selection — cross-validation (band ±1 SD, 5 folds)", fontsize=10.5, fontweight="bold", loc="left", pad=6)
    lbl(axe, "e")

    # f — parsimony, held-out test (confirmation vs competitors)
    axf = fig.add_subplot(gs[2, 1])
    axf.axvline(450, color="#ccc", lw=1.0, ls=":", zorder=1)
    axf.axhline(oai_em, color=BLUE, lw=1.3, ls=(0, (5, 4)), zorder=2, label=f"OligoAI ({oai_em:.2f})")
    axf.axhline(ow_em, color=TEAL, lw=1.3, ls=(0, (5, 4)), zorder=2, label=f"{ow_lbl} ({ow_em:.2f})")
    axf.plot(ptK, ptm, color=ACCENT, lw=2.0, zorder=4, label=f"TAUSO (K=450: {ptmap[450]:.2f})")
    droplabels(axf, lambda k: ptmap[k], ACCENT)
    logx(axf); axf.set_ylabel("within-experiment Spearman", fontsize=10)
    axf.set_title("Parsimony — held-out test (confirmation)", fontsize=10.5, fontweight="bold", loc="left", pad=6)
    axf.legend(frameon=False, fontsize=8.5, loc="center right")
    lbl(axf, "f")

fig.subplots_adjust(top=0.915 if WITH_PARSIMONY else 0.885, bottom=0.045 if WITH_PARSIMONY else 0.065,
                    left=0.085, right=0.975, hspace=0.5 if WITH_PARSIMONY else 0.42, wspace=0.32)
fig.suptitle("TAUSO vs. the field — held-out performance", fontsize=16, fontweight="bold",
             x=0.012, y=0.975 if WITH_PARSIMONY else 0.965, ha="left")
fig.text(0.012, 0.944 if WITH_PARSIMONY else 0.928, f"strict 5-10-5 MOE / 3-10-3 cEt gapmers, OligoAI-common test rows: "
         f"n = {len(S):,} ASOs / {n_cid} patents / {n_gxc} gene×cell cohorts", fontsize=9.5, color="#666")
cap = (
    "# Fig 2 — TAUSO vs. the field (held-out performance)\n\n"
    "Held-out **test** split, strict 5-10-5 2′-MOE / 3-10-3 cEt gapmers, rows OligoAI also scored "
    f"(n = {len(S):,} ASOs / {n_cid} patents / {n_gxc} gene×cell cohorts). TAUSO = the held-out model "
    "(trained on train+val, scored once), not the all-data deployment model.\n\n"
    "**(a, b)** Median within-experiment and within-cohort Spearman vs measured inhibition; TAUSO leads both. "
    "**(c)** Per experiment, take each method's top-5%-ranked ASOs and measure their actual % knockdown "
    "(mean within the experiment); the **violin is the full distribution across experiments** (thick bar = "
    "inter-quartile range, dot = median), from the random floor (a random 5% ≈ the experiment mean) to the "
    "**best possible** (the true best 5%). OligoWalk·intra (thermodynamic self-structure) is included as a baseline, "
    "shown only on the experiments where it is defined (n annotated) and sitting near the random floor. "
    "TAUSO's advantage over OligoAI is in the **floor, not the median** "
    "(70 vs 67%): its lower quartile is 61 vs 52% and it yields a failing shortlist (<50% knockdown) in only "
    "9% of experiments vs 21% (paired Wilcoxon p = 1.3×10⁻⁸; TAUSO's shortlist beats OligoAI's in 63%). "
    "**(d)** Per-chemistry within-experiment Spearman — TAUSO's margin over OligoAI is largest on cEt, whose "
    "sugar OligoAI cannot parameterise.")
if WITH_PARSIMONY:
    cap += (" **(e, f)** Feature parsimony: each K trains on the top-K RMSE-ranked "
            "feature set. (e) nested gene-grouped **cross-validation** over train+val (band ±1 SD over 5 folds) — the "
            f"selection curve, ★ = chosen K=450; (f) **held-out test** (deployed feature order, K=450 reproduces the "
            f"shipped model {ptmap[450]:.2f}) vs OligoAI ({oai_em:.2f}) and {ow_lbl} ({ow_em:.2f}). Both flat to ≈50 "
            "features; low-K drop-off values labelled.")
cap += " Draft."
save(fig, CAT, "fig2_performance_abIQR" if AB_IQR else "fig2_performance", caption=cap)
print(f"fig2 built | n={len(S)} cid_TAUSO={cid['TAUSO']:.3f} top5(median): " +
      " ".join(f"{k}={np.median(t5[k]):.0f}" for k in names) +
      (f" | parsimony test K450 {ptmap[450]:.3f}" if WITH_PARSIMONY else ""))

# ---- opt-in standalone comparison of two IQR-labelling styles (FIG2C_COMPARE=1) ----
# A: TAUSO-only median+IQR as light reference lines (current panel c).  B: every method labelled with its own median+IQR.
if os.environ.get("FIG2C_COMPARE"):
    def _violins(ax, mode):
        for i, (n, c) in enumerate(zip(names, cols)):
            lo, hi = t5[n].min(), t5[n].max(); mm = (ygrid >= lo) & (ygrid <= hi)
            w = dens[n] / kmax * 0.42
            ax.fill_betweenx(ygrid[mm], i - w[mm], i + w[mm], facecolor=c, alpha=0.5, edgecolor=c, linewidth=0.8, zorder=3)
            qa, med, qb = np.percentile(t5[n], [25, 50, 75])
            ax.plot([i, i], [qa, qb], color="#222", lw=5, solid_capstyle="butt", zorder=4)
            ax.scatter([i], [med], color="white", s=14, edgecolor="#222", lw=0.5, zorder=5)
            if mode == "all":                                   # label median + IQR for EVERY method
                ax.text(i + 0.13, med, f"{med:.1f}", color=INK_, fontsize=7.2, va="center", ha="left",
                        fontweight="bold", zorder=7, path_effects=_halo)
                for q in (qa, qb):
                    ax.text(i + 0.13, q, f"{q:.1f}", color="#9aa3b0", fontsize=5.8, va="center", ha="left",
                            zorder=7, path_effects=_halo)
        if mode == "tauso":                                     # only TAUSO's median + IQR, as light reference lines
            ti = names.index("TAUSO"); q1, md, q3 = np.percentile(t5["TAUSO"], [25, 50, 75])
            ax.axhline(md, color=ACCENT, ls="-", lw=1.0, alpha=0.40, zorder=2)
            for q in (q1, q3): ax.axhline(q, color=ACCENT, ls="--", lw=0.7, alpha=0.28, zorder=2)
            ax.text(ti + 0.13, md, f"{md:.1f}", color=INK_, fontsize=8.5, va="center", ha="left",
                    fontweight="bold", zorder=7, path_effects=_halo)
            for q in (q1, q3):
                ax.text(ti + 0.13, q, f"{q:.1f}", color="#9aa3b0", fontsize=6.5, va="center", ha="left",
                        zorder=7, path_effects=_halo)
        ax.set_xticks(xs); ax.set_xticklabels([n.replace(" ", "\n").replace("·", "\n·") for n in names], fontsize=8)
        ax.set_xlim(-0.6, len(names) - 0.4); ax.set_ylim(0, 105)
        ax.set_ylabel("top-5% knockdown (%), per experiment", fontsize=9)
    cfig, (cax1, cax2) = plt.subplots(1, 2, figsize=(13.4, 5.6))
    _violins(cax1, "tauso"); cax1.set_title("A · TAUSO median + IQR as reference lines", fontsize=10.5, fontweight="bold", loc="left", pad=6)
    _violins(cax2, "all");   cax2.set_title("B · every method labelled (median + IQR)", fontsize=10.5, fontweight="bold", loc="left", pad=6)
    cfig.suptitle("Top-5% knockdown — IQR-labelling options, side by side", fontsize=12.5, fontweight="bold", x=0.015, ha="left")
    cfig.tight_layout(rect=[0, 0, 1, 0.95])
    _cmp = Path(__file__).resolve().parents[1] / "out" / CAT / "fig2c_iqr_compare.png"
    cfig.savefig(_cmp, dpi=200, bbox_inches="tight"); print(f"wrote {_cmp}")
