"""
Signal test for two candidate sequence features, run univariately against Inhibition(%):

  1. ASO self G-quadruplex propensity (G4Hunter) -- hypothesis: flags POOR performers
     better than the existing Sequence_ggg_counts. Evaluated as poor-performer detection.
  2. Modification-weighted self-structure ("stem-in-wing"): fold the ASO, ask how much of
     the predicted stem sits in the high-affinity wings. Evaluated as median per-cohort
     (custom_id) Spearman, the ASO-selection metric.

No model is trained; these are univariate diagnostics. Run:
    python -u -m notebooks.models.seq_feature_signal_g4_selfstruct
"""
import logging
import re
import time

import numpy as np
import pandas as pd
import RNA
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

from notebooks.models.utility import load_and_validate_final_data
from tauso.data.consts import CHEMICAL_PATTERN, INHIBITION, SEQUENCE

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

COHORT = "custom_id"
MIN_COHORT = 8          # min rows for a stable per-cohort Spearman
_G4_CANON = re.compile(r"G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,}")
_G4_RELAX = re.compile(r"G{2,}\w{1,12}G{2,}\w{1,12}G{2,}\w{1,12}G{2,}")


# ----------------------------- G4Hunter -----------------------------
def _g4hunter_scores(seq):
    """Signed per-position G4Hunter scores: G-runs -> +min(run,4), C-runs -> -min(run,4)."""
    seq = seq.upper().replace("U", "T")
    n = len(seq)
    scores = np.zeros(n)
    i = 0
    while i < n:
        b = seq[i]
        if b in "GC":
            j = i
            while j < n and seq[j] == b:
                j += 1
            val = min(j - i, 4) * (1 if b == "G" else -1)
            scores[i:j] = val
            i = j
        else:
            i += 1
    return scores


def g4_features(seq):
    s = _g4hunter_scores(seq)
    n = len(s)
    if n == 0:
        return dict(g4h_mean=0.0, g4h_max5=0.0, n_gtract3=0, g4_canon=0, g4_relax=0)
    w = min(5, n)
    max5 = max(s[k:k + w].mean() for k in range(0, n - w + 1))
    up = seq.upper().replace("U", "T")
    return dict(
        g4h_mean=float(s.mean()),
        g4h_max5=float(max5),
        n_gtract3=len(re.findall(r"G{3,}", up)),
        g4_canon=int(bool(_G4_CANON.search(up))),
        g4_relax=int(bool(_G4_RELAX.search(up))),
    )


# ------------------------- stem-in-wing self-structure -------------------------
def selfstruct_features(seq, pattern):
    """Fold the ASO; measure how much of the predicted stem sits in the modified wings.
    Wing = non-'d' position in the chemical pattern (high-affinity sugar); gap = 'd' (DNA)."""
    struct, mfe = RNA.fold(seq.upper().replace("U", "T"))
    paired = np.array([c != "." for c in struct])
    n_pair = int(paired.sum())
    out = dict(internal_mfe=float(mfe), paired_frac=n_pair / len(struct) if struct else 0.0)
    if pattern is None or not isinstance(pattern, str) or len(pattern) != len(struct) or n_pair == 0:
        out["stem_wing_frac"] = 0.0
        out["selfstruct_wing"] = 0.0
        return out
    wing = np.array([c != "d" for c in pattern.lower()])
    wing_frac = float((paired & wing).sum()) / n_pair
    out["stem_wing_frac"] = wing_frac
    out["selfstruct_wing"] = abs(mfe) * wing_frac      # stable self-structure parked in sticky chemistry
    return out


# ----------------------------- evaluation -----------------------------
def median_cohort_spearman(df, col, target=INHIBITION):
    rhos = []
    for _, g in df.groupby(COHORT):
        if g[col].nunique() < 3 or len(g) < MIN_COHORT:
            continue
        rho, _ = spearmanr(g[col], g[target])
        if np.isfinite(rho):
            rhos.append(rho)
    return (np.median(rhos), np.mean(rhos), len(rhos))


def poor_performer_report(df, cols, poor_q=0.25):
    """Higher predictor should mean worse ASO. ROC-AUC of predictor vs is-poor (bottom quartile)."""
    thr = df[INHIBITION].quantile(poor_q)
    is_poor = (df[INHIBITION] <= thr).astype(int)
    base_rate = is_poor.mean()
    logger.info("\n=== POOR-PERFORMER DETECTION (poor = bottom %.0f%% Inhibition, <=%.1f; base rate %.3f, n=%d) ===",
                poor_q * 100, thr, base_rate, len(df))
    logger.info("%-16s %8s %8s %10s %10s %12s", "feature", "AUC", "globRho", "n(>0)", "mean_inh>0", "lift_poor>0")
    for c in cols:
        x = df[c].astype(float)
        auc = roc_auc_score(is_poor, x) if x.nunique() > 1 else float("nan")
        rho, _ = spearmanr(x, df[INHIBITION])
        pos = x > x.min()
        n_pos = int(pos.sum())
        mean_inh_pos = df.loc[pos, INHIBITION].mean() if n_pos else float("nan")
        lift = (is_poor[pos].mean() / base_rate) if n_pos else float("nan")
        logger.info("%-16s %8.3f %8.3f %10d %10.1f %12.2f", c, auc, rho, n_pos, mean_inh_pos, lift)


def main():
    final_data, _ = load_and_validate_final_data(version="oligo")
    df = final_data.dropna(subset=[INHIBITION, SEQUENCE]).copy()
    logger.info("Loaded %d rows, %d cohorts (%s)", len(df), df[COHORT].nunique(), COHORT)

    # Compute on unique sequences / (sequence, pattern) pairs, then map back.
    t0 = time.time()
    uniq_seq = pd.Index(df[SEQUENCE].unique())
    g4 = pd.DataFrame([g4_features(s) for s in uniq_seq], index=uniq_seq)
    for c in g4.columns:
        df[c] = df[SEQUENCE].map(g4[c])

    pairs = df[[SEQUENCE, CHEMICAL_PATTERN]].drop_duplicates()
    ss = pd.DataFrame([selfstruct_features(r[SEQUENCE], r[CHEMICAL_PATTERN]) for _, r in pairs.iterrows()])
    ss.index = pd.MultiIndex.from_frame(pairs)
    key = pd.MultiIndex.from_frame(df[[SEQUENCE, CHEMICAL_PATTERN]])
    for c in ss.columns:
        df[c] = ss[c].reindex(key).to_numpy()
    logger.info("Feature compute: %.1fs (%d uniq seqs, %d uniq seq/pattern)", time.time() - t0, len(uniq_seq), len(pairs))

    # ---- (1) G4 vs poor performers, head-to-head with the existing ggg feature ----
    g4_cols = ["g4h_mean", "g4h_max5", "n_gtract3", "g4_canon", "g4_relax"]
    if "Sequence_ggg_counts" in df.columns:
        g4_cols.append("Sequence_ggg_counts")
    poor_performer_report(df, g4_cols, poor_q=0.25)
    poor_performer_report(df, g4_cols, poor_q=0.10)

    # ---- (2) stem-in-wing self-structure: median per-cohort Spearman ----
    logger.info("\n=== MEDIAN PER-COHORT SPEARMAN vs Inhibition (cohorts >= %d rows) ===", MIN_COHORT)
    logger.info("%-18s %10s %10s %8s", "feature", "median_rho", "mean_rho", "n_coh")
    ss_cols = ["selfstruct_wing", "stem_wing_frac", "internal_mfe", "paired_frac"]
    if "Sequence_internal_fold" in df.columns:
        ss_cols.append("Sequence_internal_fold")
    for c in ss_cols:
        med, mean, ncoh = median_cohort_spearman(df, c)
        logger.info("%-18s %10.4f %10.4f %8d", c, med, mean, ncoh)


if __name__ == "__main__":
    main()
