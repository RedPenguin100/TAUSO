"""Shared helpers for ASO model training and evaluation: data splitting + ranking metrics."""
import numpy as np
from scipy.stats import spearmanr


def split_data(df, features, split_col="split", train_val="train", val_val="val", test_val="test"):
    return (
        df[df[split_col] == train_val].copy(),
        df[df[split_col] == val_val].copy(),
        df[df[split_col] == test_val].copy(),
    )


def get_large_cohort_indices(df, group_cols, min_size=50):
    df = df.reset_index(drop=True)
    return [g.index.values for _, g in df.groupby(group_cols) if len(g) >= min_size]


def _to_numpy(arr):
    """Convert to NumPy, handling CuPy arrays returned by GPU prediction."""
    if hasattr(arr, "get"):
        return arr.get()
    return np.asarray(arr)


def _collect_group_stats(y_true, preds, eval_groups):
    """Per-group (spearman list, top-1% arrays, top-5% arrays)."""
    preds = _to_numpy(preds)
    y_true = _to_numpy(y_true)
    spearmans, top1s, top5s = [], [], []
    for idxs in eval_groups:
        t, p = y_true[idxs], preds[idxs]
        corr, _ = spearmanr(t, p)
        if not np.isnan(corr):
            spearmans.append(corr)
        n = len(t)
        k1, k5 = max(1, int(n * 0.01)), max(1, int(n * 0.05))
        top1s.append(t[np.argpartition(p, -k1)[-k1:]])
        top5s.append(t[np.argpartition(p, -k5)[-k5:]])
    return spearmans, top1s, top5s


def calculate_metrics(preds, y_true, eval_groups):
    spearmans, top1s, top5s = _collect_group_stats(y_true, preds, eval_groups)
    return {
        "spearman": float(np.nanmedian(spearmans)) if spearmans else 0.0,
        "top1_inhibition": float(np.nanmean([np.mean(t) for t in top1s])) if top1s else 0.0,
        "top5_inhibition": float(np.nanmean([np.mean(t) for t in top5s])) if top5s else 0.0,
    }
