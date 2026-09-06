"""Train the hl-LOW/rmse config on the current feature set.

  python notebooks/models/train_hl_low_rmse.py              # train+val, seed 1, score the test split
  python notebooks/models/train_hl_low_rmse.py --data all   # train on train+val+test

Writes the booster and the feature list it was trained on next to the other model artefacts.
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # repo root, for the notebooks.* imports
from notebooks.models import common

CONFIG = json.loads((Path(__file__).parent / "hl_low_rmse_parameters.json").read_text())["hl_low_rmse"]
LABEL = "hl_low_rmse"


def current_pipeline_features():
    """Names the feature pipeline writes a shard for, i.e. what it still computes."""
    from notebooks.features.feature_extraction import _get_saved_features_dir
    from tauso.populate.feature_cache import loose_shard_dir

    shard_dir = Path(loose_shard_dir(_get_saved_features_dir("oligo")))
    return {p.stem for p in shard_dir.iterdir() if p.suffix in (".parquet", ".csv")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", choices=["trainval", "all"], default="trainval")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    df, features = common.load_dataset()
    print(f"loaded {len(df)} rows, {len(features)} features", flush=True)

    # The loader fills anything the per-feature shards do not cover from the wide parquet cache,
    # which still holds columns the pipeline has since stopped producing. Training on those would
    # make the model unscorable, so the shards are what defines the feature set.
    produced = current_pipeline_features()
    retired = sorted(set(features) - produced)
    features = [f for f in features if f in produced]
    print(f"{len(retired)} features come only from the wide cache and are dropped: {retired}", flush=True)

    features = [f for f in features if df[f].nunique(dropna=True) > 1]  # a constant column carries no split
    trv, test = common.split(df)
    train_df = trv if args.data == "trainval" else df

    params, rounds, variant = CONFIG["params"], CONFIG["num_boost_round"], CONFIG["variant"]
    print(f"training {LABEL} ({variant}) on {args.data}: {len(train_df)} rows, "
          f"{len(features)} features, {rounds} rounds, seed {args.seed}", flush=True)

    model = common.train(train_df, features, variant, params, rounds, seed=args.seed)

    stem = common.RESULTS_DIR / f"model_{LABEL}_{args.data}_seed{args.seed}"
    stem.parent.mkdir(parents=True, exist_ok=True)
    model.save_model(str(stem.with_suffix(".json")))
    stem.with_suffix(".features.txt").write_text("\n".join(features) + "\n")

    if args.data == "trainval":
        scores = common.metrics_on(model, test, features)
        print("\nTEST (held out, n=%d)" % len(test))
        for metric in common.METRICS:
            print(f"  {metric:>8}: {scores[metric]:.4f}")

    print(f"\n  -> {stem.with_suffix('.json')}\n  -> {stem.with_suffix('.features.txt')}")


if __name__ == "__main__":
    main()
