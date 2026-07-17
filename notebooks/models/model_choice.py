"""Model choice -- CV reaffirmation.

Cross-validates the v15 HPO champions (best_models_parameters.json) on train+val with the frozen
gene-grouped folds and reports the full metric suite per model, so we can see which config is best.
Each model trains with its own objective (clean_exp / clean_gene / reg) and tuned params.

CV only -- the held-out test is deliberately untouched here (test validation / deploy is deploy.py).

Run:  python notebooks/models/model_choice.py
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # repo root, for the notebooks.* imports
from notebooks.models import common

BEST_MODELS = json.loads((Path(__file__).parent / "best_models_parameters.json").read_text())


def cv_scores(trv, features, variant, params, rounds):
    """Mean +/- std of the full metric suite over the frozen gene-grouped CV folds."""
    return common.mean_std([common.metrics_on(common.train(tr, features, variant, params, rounds), va, features)
                            for tr, va in common.gene_cv_folds(trv)])


def main():
    df, features = common.load_dataset()
    features = [f for f in features if df[f].nunique(dropna=True) > 1]   # drop zero-variance -> the deployed 485-feature set
    trv, _ = common.split(df)                                   # train+val ONLY; test is never used for CV
    print(f"model choice (CV reaffirmation): {len(features)} feats, {len(trv)} train+val rows, "
          f"{len(BEST_MODELS)} models x 5 gene-folds", flush=True)

    cv = {}
    for k, (name, spec) in enumerate(BEST_MODELS.items(), 1):
        cv[name] = cv_scores(trv, features, spec["variant"], spec["params"], spec["num_boost_round"])
        print(f"[{k}/{len(BEST_MODELS)}] {name} ({spec['variant']}) done -- "
              f"CV exp_med {cv[name]['exp_med'][0]:.4f}  gxc_med {cv[name]['gxc_med'][0]:.4f}", flush=True)

    lines = [f"Model choice -- CV reaffirmation of the {len(BEST_MODELS)} v15 HPO champions   {len(features)} features"
             f"   train+val={len(trv)}   cv=5-fold gene-grouped   (NO held-out test)", ""]
    lines += common.metric_table("CROSS-VALIDATION (mean +/- std over 5 gene-folds)", cv) + [""]
    lines += ["best by CV metric:"]
    for m in ("exp_med", "gxc_med"):
        best = max(BEST_MODELS, key=lambda n: cv[n][m][0])
        lines.append(f"  {m:8s}: {best}  ({cv[best][m][0]:.4f} ± {cv[best][m][1]:.4f})")
    lines += ["", "models (objective; sweep CV for provenance):"]
    for name, spec in BEST_MODELS.items():
        lines.append(f"  {name:20s} {spec['variant']:11s}  sweep CV exp_med {spec['cv_exp_med']:.4f}  gxc_med {spec['cv_gxc_med']:.4f}")
    common.write_report(lines, "model_choice.txt")

    # Full result, machine-readable: every model x metric, as (mean, std) over the CV folds.
    result = {
        "n_features": len(features), "train_val_rows": len(trv), "cv_folds": 5,
        "metrics": list(common.METRICS),
        "cv": {name: {m: list(cv[name][m]) for m in common.METRICS} for name in BEST_MODELS},
        "models": {name: {"variant": spec["variant"], "params": spec["params"],
                          "num_boost_round": spec["num_boost_round"]} for name, spec in BEST_MODELS.items()},
    }
    (common.RESULTS_DIR / "model_choice.json").write_text(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
