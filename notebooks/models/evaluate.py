"""Hidden-test metrics, all against actual inhibition: within-experiment (custom_id) and cohort (gene x cell)
Spearman, top-5% predicted mean actual inhibition, precision@k, global Pearson, RMSE/MAE."""
import numpy as np
from scipy.stats import pearsonr, spearmanr


def _grouped_spearman(pred, y, group, min_size=5):
    out = []
    for g in np.unique(group):
        m = group == g
        if m.sum() >= min_size and len(np.unique(y[m])) > 1:
            r = spearmanr(pred[m], y[m]).correlation
            if np.isfinite(r):
                out.append(r)
    return out


def evaluate(pred, y, custom_id, cohort):
    exp = _grouped_spearman(pred, y, custom_id)
    gxc = _grouped_spearman(pred, y, cohort)
    top5, p5, p10 = [], [], []
    for c in np.unique(custom_id):
        m = custom_id == c
        n = m.sum()
        if n >= 20 and len(np.unique(y[m])) > 1:
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
