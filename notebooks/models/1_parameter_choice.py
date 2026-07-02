"""Step 1 -- parameter choice: default vs regularized XGBoost, by gene-grouped CV on train+val."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # repo root, for the notebooks.* imports
from notebooks.models import common

# default = the stock XGBoost values made explicit (so the report states them); custom = our tuned config.
DEFAULT = dict(tree_method="hist", device="cuda", max_depth=6, learning_rate=0.3,
               subsample=1.0, colsample_bytree=1.0, min_child_weight=1, reg_lambda=1.0, reg_alpha=0.0)
CONFIGS = {"default": (DEFAULT, 100), "custom": (common.PARAMS, common.ROUNDS)}   # (params, num_boost_round)
HEADLINE = "reg"
OBJECTIVES = ["reg", "clean_exp", "clean_gene", "pair_exp", "pair_gene"]


def cv_scores(trv, features, variant, params, rounds):
    """Mean +/- std of the metric suite over the frozen gene-grouped CV folds."""
    return common.mean_std([common.metrics_on(common.train(tr, features, variant, params, rounds), va, features)
                            for tr, va in common.gene_cv_folds(trv)])


def _config_line(name):
    params, rounds = CONFIGS[name]
    return f"  {name:8s}: " + ", ".join(f"{k}={v}" for k, v in params.items()) + f", num_boost_round={rounds}"


def main():
    df, features = common.load_dataset()
    trv, _ = common.split(df)

    cv = {v: {name: cv_scores(trv, features, v, params, rounds) for name, (params, rounds) in CONFIGS.items()}
          for v in OBJECTIVES}

    lines = [f"Parameter choice -- default vs custom XGBoost   {len(features)} features"
             f"   train+val={len(trv)}   cv=5-fold gene-grouped", ""]
    lines += common.metric_table(f"[{HEADLINE}] CROSS-VALIDATION (mean +/- std over 5 gene-folds)", cv[HEADLINE]) + [""]
    lines += ["ALL OBJECTIVES -- CV exp_med:",
              f"  {'objective':12s} {'default':>12} {'custom':>12} {'delta':>9}"]
    for v in OBJECTIVES:
        d, c = cv[v]["default"]["exp_med"], cv[v]["custom"]["exp_med"]
        lines.append(f"  {v:12s} {d[0]:6.3f}±{d[1]:.3f} {c[0]:6.3f}±{c[1]:.3f} {c[0] - d[0]:+9.3f}")
    lines += ["", "configs:", _config_line("default"), _config_line("custom")]
    common.write_report(lines, "parameter_choice.txt")


if __name__ == "__main__":
    main()
