"""Do the literature-grounded toxicity/liability motifs predict POOR-performing ASOs?

Univariate + composite diagnostics against Inhibition(%). Reminder: toxicity != inefficacy,
so hepatotox (TCC/TGC) and TLR9/CpG motifs are NOT expected to predict knockdown; only the
G-quadruplex / G-run family has a mechanistic route to reduced efficacy. We therefore also
report mean inhibition among flagged oligos to show whether "liable" ASOs are still effective.

    python -u -m notebooks.models.toxicity_motif_signal
"""
import logging

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

from tauso.data.consts import CHEMICAL_PATTERN, INHIBITION, SEQUENCE
from tauso.features.sequence.seq_features import count_g_runs
from tauso.features.sequence.toxicity_features import (
    cpg_count,
    g4hunter_max,
    g4hunter_mean,
    hepatotox_gap_motif_count,
    hepatotox_motif_count,
    max_g_run,
    three_prime_g_fraction,
    three_prime_g_free_length,
    tlr9_human_motif_count,
)

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)
COHORT = "custom_id"
MIN_COHORT = 8
DATA = "notebooks/data/OligoAI/raw_data/aso_inhibitions_with_canonical_gene_processed_averaged.csv.gz"

# group -> (feature name, function-of-row)
TOX = {  # pure-toxicity axis (not expected to predict efficacy)
    "hepatotox_gap": lambda s, p: hepatotox_gap_motif_count(s, p),
    "hepatotox_whole": lambda s, p: hepatotox_motif_count(s),
    "cpg_count": lambda s, p: cpg_count(s),
    "tlr9_GTCGTT": lambda s, p: tlr9_human_motif_count(s),
    "tp_g_frac": lambda s, p: three_prime_g_fraction(s),
    "tp_g_free_len": lambda s, p: three_prime_g_free_length(s),
}
GFAM = {  # G-quadruplex / aggregation axis (the efficacy-relevant one)
    "g4hunter_mean": lambda s, p: g4hunter_mean(s),
    "g4hunter_max": lambda s, p: g4hunter_max(s),
    "max_g_run": lambda s, p: max_g_run(s),
    "ggg_counts(existing)": lambda s, p: count_g_runs(s),
}


def _z(s):
    sd = s.std()
    return (s - s.mean()) / sd if sd else s * 0.0


def poor_report(df, cols, q=0.25):
    thr = df[INHIBITION].quantile(q)
    is_poor = (df[INHIBITION] <= thr).astype(int)
    base = is_poor.mean()
    logger.info("\n=== POOR-PERFORMER DETECTION (poor = bottom %.0f%% Inhibition <= %.1f; base %.3f) ===",
                q * 100, thr, base)
    logger.info("%-22s %7s %8s %9s %12s %11s", "feature", "AUC", "globRho", "n(>0)", "meanInh(>0)", "lift_poor")
    for c in cols:
        x = df[c].astype(float)
        auc = roc_auc_score(is_poor, x) if x.nunique() > 1 else float("nan")
        rho = spearmanr(x, df[INHIBITION]).statistic
        pos = x > x.min()
        npos = int(pos.sum())
        mean_pos = df.loc[pos, INHIBITION].mean() if npos else float("nan")
        lift = (is_poor[pos].mean() / base) if npos else float("nan")
        logger.info("%-22s %7.3f %8.3f %9d %12.1f %11.2f", c, auc, rho, npos, mean_pos, lift)


def med_cohort_spearman(df, cols):
    logger.info("\n=== MEDIAN PER-COHORT SPEARMAN vs Inhibition (cohorts >= %d rows) ===", MIN_COHORT)
    logger.info("%-22s %10s %8s", "feature", "median_rho", "n_coh")
    for c in cols:
        rhos = []
        for _, g in df.groupby(COHORT):
            if len(g) < MIN_COHORT or g[c].nunique() < 3:
                continue
            r = spearmanr(g[c], g[INHIBITION]).statistic
            if np.isfinite(r):
                rhos.append(r)
        logger.info("%-22s %10.4f %8d", c, np.median(rhos) if rhos else float("nan"), len(rhos))


def main():
    df = pd.read_csv(DATA, usecols=[SEQUENCE, CHEMICAL_PATTERN, INHIBITION, COHORT]).dropna(subset=[INHIBITION, SEQUENCE])
    logger.info("Loaded %d rows, %d cohorts. Overall mean Inhibition = %.1f", len(df), df[COHORT].nunique(), df[INHIBITION].mean())

    pairs = df[[SEQUENCE, CHEMICAL_PATTERN]].drop_duplicates()
    allspec = {**TOX, **GFAM}
    feats = pd.DataFrame({name: [fn(r[SEQUENCE], r[CHEMICAL_PATTERN]) for _, r in pairs.iterrows()]
                          for name, fn in allspec.items()})
    feats.index = pd.MultiIndex.from_frame(pairs)
    key = pd.MultiIndex.from_frame(df[[SEQUENCE, CHEMICAL_PATTERN]])
    for c in feats.columns:
        df[c] = feats[c].reindex(key).to_numpy()

    # composites (z-scored sums)
    df["TOX_sum_z"] = sum(_z(df[c]) for c in TOX)
    df["GFAM_sum_z"] = sum(_z(df[c]) for c in GFAM)
    df["ALL_sum_z"] = df["TOX_sum_z"] + df["GFAM_sum_z"]

    order = list(TOX) + list(GFAM) + ["TOX_sum_z", "GFAM_sum_z", "ALL_sum_z"]
    poor_report(df, order, q=0.25)
    poor_report(df, order, q=0.10)
    med_cohort_spearman(df, order)

    # tox != inefficacy: are flagged-as-liable oligos still effective?
    logger.info("\n=== TOX != INEFFICACY CHECK (mean Inhibition: flagged vs not) ===")
    for c in ["hepatotox_gap", "cpg_count", "tlr9_GTCGTT", "max_g_run"]:
        flag = df[c] > (2 if c == "max_g_run" else 0)
        logger.info("  %-16s flagged n=%6d meanInh=%.1f | not-flagged meanInh=%.1f",
                    c, int(flag.sum()), df.loc[flag, INHIBITION].mean(), df.loc[~flag, INHIBITION].mean())


if __name__ == "__main__":
    main()
