"""Script 1 -- model choice (CV reaffirmation).

Cross-validates the six v15 HPO champions (champions.CHAMPIONS) on train+val with the frozen
gene-grouped folds and reports the full metric suite per champion, so we can see which config is
best. Each champion trains with its own objective (clean_exp / clean_gene / reg) and tuned params.

CV only -- the held-out test is deliberately untouched here (test validation / deploy is deploy.py).

Run:  python notebooks/models/1_model_choice.py
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # repo root, for the notebooks.* imports
from notebooks.models import champions, common


def cv_scores(trv, features, variant, params, rounds):
    """Mean +/- std of the full metric suite over the frozen gene-grouped CV folds."""
    return common.mean_std([common.metrics_on(common.train(tr, features, variant, params, rounds), va, features)
                            for tr, va in common.gene_cv_folds(trv)])


def main():
    df, features = common.load_dataset()
    features = [f for f in features if df[f].nunique(dropna=True) > 1]   # drop zero-variance -> the deployed 485-feature set
    trv, _ = common.split(df)                                   # train+val ONLY; test is never used for CV
    names = list(champions.CHAMPIONS)
    print(f"script 1 (CV reaffirmation): {len(features)} feats, {len(trv)} train+val rows, "
          f"{len(names)} champions x 5 gene-folds", flush=True)

    cv = {}
    for k, name in enumerate(names, 1):
        params, rounds, variant = champions.config(name)
        cv[name] = cv_scores(trv, features, variant, params, rounds)
        print(f"[{k}/{len(names)}] {name} ({variant}) done -- "
              f"CV exp_med {cv[name]['exp_med'][0]:.4f}  gxc_med {cv[name]['gxc_med'][0]:.4f}", flush=True)

    lines = [f"Model choice -- CV reaffirmation of the {len(names)} v15 HPO champions   {len(features)} features"
             f"   train+val={len(trv)}   cv=5-fold gene-grouped   (NO held-out test)", ""]
    lines += common.metric_table("CROSS-VALIDATION (mean +/- std over 5 gene-folds)", cv) + [""]
    lines += ["best by CV metric:"]
    for m in ("exp_med", "gxc_med"):
        best = max(names, key=lambda n: cv[n][m][0])
        lines.append(f"  {m:8s}: {best}  ({cv[best][m][0]:.4f} ± {cv[best][m][1]:.4f})")
    lines += ["", "champions (objective; sweep CV for provenance):"]
    for name in names:
        c = champions.CHAMPIONS[name]
        lines.append(f"  {name:20s} {c['variant']:11s}  sweep CV exp_med {c['cv_exp_med']:.4f}  gxc_med {c['cv_gxc_med']:.4f}")
    common.write_report(lines, "model_choice.txt")

    # Full result, machine-readable: every champion x metric, as (mean, std) over the CV folds.
    result = {
        "n_features": len(features), "train_val_rows": len(trv), "cv_folds": 5,
        "metrics": list(common.METRICS),
        "cv": {name: {m: list(cv[name][m]) for m in common.METRICS} for name in names},
        "champions": {name: {"variant": champions.CHAMPIONS[name]["variant"],
                             "params": champions.CHAMPIONS[name]["params"],
                             "num_boost_round": champions.CHAMPIONS[name]["rounds"]} for name in names},
    }
    (common.RESULTS_DIR / "model_choice.json").write_text(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
