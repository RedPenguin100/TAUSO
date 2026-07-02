"""Hidden-test metrics, all against actual inhibition: within-experiment (custom_id) and cohort (gene x cell)
Spearman, top-5% predicted mean actual inhibition, precision@k, global Pearson, RMSE/MAE."""
import os

import numpy as np
from scipy.stats import pearsonr, spearmanr

# Minimum ASOs in a group to score it. Set once per whole run via the environment so every step agrees.
MIN_GROUP_N = int(os.environ.get("TAUSO_MIN_GROUP_N", 5))              # per-group Spearman (custom_id / cohort_id)
MIN_GROUP_N_TOPK = int(os.environ.get("TAUSO_MIN_GROUP_N_TOPK", 20))  # top-5% mean inhib / precision@k (custom_id)


def grouped_spearman(pred, y, groups, min_n=MIN_GROUP_N, return_sizes=False):
    """Spearman(pred, y) within each group of >= min_n members and >1 distinct label. A group whose
    PREDICTIONS are all tied (undefined correlation) scores 0 -- a failed ranking is penalised, not
    dropped, so every model/criterion is scored on the same data-determined set of groups. Groups whose
    TRUE labels are all tied are unrankable and excluded. Returns correlations (+ sizes if requested)."""
    corrs, sizes = [], []
    for g in np.unique(groups):
        m = groups == g
        n = int(m.sum())
        if n >= min_n and len(np.unique(y[m])) > 1:
            r = spearmanr(pred[m], y[m]).correlation if len(np.unique(pred[m])) > 1 else 0.0
            corrs.append(r if np.isfinite(r) else 0.0)
            sizes.append(n)
    return (np.array(corrs), np.array(sizes)) if return_sizes else np.array(corrs)


def evaluate(pred, y, custom_id, cohort):
    exp = grouped_spearman(pred, y, custom_id)
    gxc = grouped_spearman(pred, y, cohort)
    top5, p5, p10 = [], [], []
    for c in np.unique(custom_id):
        m = custom_id == c
        n = m.sum()
        if n >= MIN_GROUP_N_TOPK and len(np.unique(y[m])) > 1:
            order, truth = np.argsort(-pred[m]), np.argsort(-y[m])
            k = max(1, int(np.ceil(0.05 * n)))
            top5.append(y[m][order[:k]].mean())
            for kk, bucket in ((5, p5), (10, p10)):
                if n >= 2 * kk:
                    bucket.append(len(set(order[:kk]) & set(truth[:kk])) / kk)
    err = pred - y
    return dict(exp_med=np.median(exp), exp_mean=np.mean(exp),
                gxc_med=np.median(gxc), gxc_mean=np.mean(gxc),
                top5=np.median(top5), p5=np.mean(p5), p10=np.mean(p10),
                globP=pearsonr(pred, y)[0],
                RMSE=np.sqrt((err ** 2).mean()), MAE=np.abs(err).mean())
