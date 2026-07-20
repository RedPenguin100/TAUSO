"""Fig 2 data preparation — all numbers, NO matplotlib.

Assembles the shared held-out comparison subset (OligoAI split, strict gapmers, rows OligoAI also
scored) and every per-panel quantity. The heavy sweeps live in models/ and are read here from their
CSVs (they are logic, run once): parsimony_gain.py -> parsimony curve, screening_effort.py -> effort.
Panels consume a single Fig2Data instance; they do no data work.
"""
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

GRAPHS = Path(__file__).resolve().parents[2]          # notebooks/graphs
MODELS = GRAPHS / "models"
sys.path.insert(0, str(MODELS))
import _modeldata as M                                # noqa: E402  (sets up tauso/notebooks paths)
from tauso.data.consts import CHEMICAL_PATTERN        # noqa: E402

MIN_ROWS = 3                                          # min ASOs/screen for a within-screen Spearman
IDX, INH, CID = M.IDX, M.INH, M.CID
PARSIMONY_CSV = MODELS / "parsimony_gain_iterative_test.csv"
EFFORT_CSV = MODELS / "screening_effort.csv"
EFFORT_K = 5


def _chem_of(p):
    p = str(p)
    if "M" in p and "C" in p: return "mixmer"
    if "M" in p: return "2'-MOE"
    if "C" in p: return "cEt"
    if p and set(p) <= set("d"): return "PS-DNA"
    return "other"


def _effort_summary(effort, methods=("TAUSO", "OligoAI", "random"), n_boot=4000):
    """Per method: (geometric-mean fold, 95% bootstrap-CI low, high), resampling screens."""
    rng = np.random.default_rng(0)
    out = {}
    for m in methods:
        x = effort[m].to_numpy(float)
        n = len(x)
        boots = np.exp(np.log(x[rng.integers(0, n, size=(n_boot, n))]).mean(axis=1))
        lo, hi = np.percentile(boots, [2.5, 97.5])
        out[m] = (float(np.exp(np.log(x).mean())), float(lo), float(hi))
    return out


def _top5_dist(S, score, s=1.0):
    """Per experiment: mean actual knockdown of the top-5%-ranked ASOs by `score` (or random/oracle)."""
    out = []
    for _, g in S.groupby(CID):
        n = len(g)
        if n >= MIN_ROWS and g[INH].nunique() > 1:
            k = max(1, int(round(0.05 * n)))
            if score == "_random":
                out.append(g[INH].mean())
            elif score == "_oracle":
                out.append(g.nlargest(k, INH)[INH].mean())
            else:
                sel = (g[score] * s).dropna().nlargest(k)
                if len(sel):
                    out.append(g.loc[sel.index, INH].mean())
    return np.array(out)


@dataclass
class Fig2Data:
    present: dict            # method label -> score column
    sign: dict               # method label -> +1/-1 orientation
    cid: dict                # (a) per-patent median Spearman
    gxc: dict                # (b) per gene x cell median Spearman
    t5: dict                 # (c) method -> array of per-experiment top-5% knockdown
    chem_res: dict           # (d) chemistry -> {TAUSO, OligoAI}
    parsimony: pd.DataFrame  # (e) K, exp_med
    oai_em: float            # (e) OligoAI reference on this subset
    ow_em: float             # (e) best rule-based reference
    ow_lbl: str
    effort: pd.DataFrame     # (f) per-screen fold: TAUSO / OligoAI / random
    effort_stat: dict        # (f) method -> (geomean fold, 95% boot-CI lo, hi) over screens
    effort_k: int
    n_aso: int
    n_cid: int
    n_gxc: int


def load():
    avg = pd.read_csv(M.AVG, usecols=[IDX, INH, CID, "split", CHEMICAL_PATTERN, M.GENE, M.CELL])
    avg["chem"] = avg[CHEMICAL_PATTERN].map(_chem_of)
    avg["is_strict"] = avg[CHEMICAL_PATTERN].isin(M.STRICT)
    df = avg
    present = {}
    for lbl, stem in M.COMP.items():
        try:
            df = df.merge(M._shard(stem), on=IDX, how="left"); present[lbl] = stem
        except FileNotFoundError:
            pass
    df = df.merge(pd.read_parquet(M.TAUSO_PRED), on=IDX, how="left"); present["TAUSO"] = "tauso"
    df = df[(df["split"] == "test") & df[INH].notna()].reset_index(drop=True)
    df["gxc"] = df[M.GENE].astype(str) + "|" + df[M.CELL].astype(str)
    S = df[df["is_strict"] & df["oligo_ai_score"].notna()].reset_index(drop=True)
    sign = M.orient(S, present)

    cid = {lbl: M.med(S, col, CID) * sign[lbl] for lbl, col in present.items()}
    gxc = {lbl: M.med(S, col, "gxc") * sign[lbl] for lbl, col in present.items()}
    t5 = {n: _top5_dist(S, sc, sg) for n, sc, sg in [
        ("Random", "_random", 1.0),
        ("OligoWalk·intra", present.get("OligoWalk·intra"), sign.get("OligoWalk·intra", 1.0)),
        ("OligoAI", "oligo_ai_score", sign["OligoAI"]),
        ("TAUSO", "tauso", sign["TAUSO"]),
        ("best possible", "_oracle", 1.0)] if sc is not None}
    chem_res = {c: {"TAUSO": M.med(S[S.chem == c], "tauso", CID),
                    "OligoAI": M.med(S[S.chem == c], "oligo_ai_score", CID)} for c in ["2'-MOE", "cEt"]}

    parsimony = pd.read_csv(PARSIMONY_CSV).astype({"K": int}).sort_values("K")
    oai_em = M.med(S, "oligo_ai_score", CID) * sign["OligoAI"]
    _ow = {lbl: M.med(S, present[lbl], CID) * sign[lbl] for lbl in present if "OligoWalk" in lbl or lbl == "sfold"}
    ow_lbl = max(_ow, key=_ow.get)
    effort = pd.read_csv(EFFORT_CSV)

    return Fig2Data(present=present, sign=sign, cid=cid, gxc=gxc, t5=t5, chem_res=chem_res,
                    parsimony=parsimony, oai_em=oai_em, ow_em=_ow[ow_lbl], ow_lbl=ow_lbl,
                    effort=effort, effort_stat=_effort_summary(effort), effort_k=EFFORT_K,
                    n_aso=len(S), n_cid=S[CID].nunique(), n_gxc=S["gxc"].nunique())
