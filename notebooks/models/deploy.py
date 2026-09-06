"""Deploy -- train and save the final shipped model.

Trains the config in deploy_parameters.json on train+val and writes the booster, the feature
list the package scores with, and the list of those features that are never missing.

  python notebooks/models/deploy.py                 # train on train+val, seed 1  (the final model)
  python notebooks/models/deploy.py --data all      # train on all data (train+val+test)

The booster is ~100 MB and stays out of git; copy it to <data_dir>/models/ to score with it.
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # repo root, for the notebooks.* imports
from notebooks.models import common
from tauso.inference.predict import DEFAULT_VERSION, MODEL_DIR, MODEL_FILES

CONFIGS = json.loads((Path(__file__).parent / "deploy_parameters.json").read_text())


def current_pipeline_features():
    """Names the feature pipeline writes a shard for, i.e. what it still computes."""
    from notebooks.features.feature_extraction import _get_saved_features_dir
    from tauso.populate.feature_cache import loose_shard_dir

    shard_dir = Path(loose_shard_dir(_get_saved_features_dir("oligo")))
    return {p.stem for p in shard_dir.iterdir() if p.suffix in (".parquet", ".csv")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", choices=["trainval", "all"], default="trainval",
                    help="train on train+val (default, the shipped model) or all data")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--version", default=DEFAULT_VERSION, help="model version to write")
    args = ap.parse_args()

    spec = CONFIGS[args.version]
    df, features = common.load_dataset()

    # The loader fills what the per-feature shards do not cover from the wide parquet cache, which
    # still holds columns the pipeline has stopped producing. A model trained on one of those cannot
    # be scored by a feature run, so the shards are what defines the feature set.
    produced = current_pipeline_features()
    stale = sorted(set(features) - produced)
    features = [f for f in features if f in produced]
    if stale:
        print(f"{len(stale)} features come only from the wide cache and are dropped: {stale}", flush=True)

    features = [f for f in features if df[f].nunique(dropna=True) > 1]   # a constant column carries no split
    trv, test = common.split(df)
    train_df = trv if args.data == "trainval" else df

    params, rounds, variant = spec["params"], spec["num_boost_round"], spec["variant"]
    print(f"deploy {args.version} ({variant}): train on {args.data} ({len(train_df)} rows), "
          f"seed {args.seed}, {rounds} rounds, {len(features)} feats", flush=True)

    model = common.train(train_df, features, variant, params, rounds, seed=args.seed)

    booster = common.RESULTS_DIR / MODEL_FILES[args.version]["filename"]
    booster.parent.mkdir(parents=True, exist_ok=True)
    model.save_model(str(booster))

    # What the package scores with: the feature order, and the features a NaN is a failure in.
    finite = sorted(f for f in features if not train_df[f].isna().any())
    (MODEL_DIR / f"tauso_score_{args.version}.features.txt").write_text("\n".join(features) + "\n")
    (MODEL_DIR / f"tauso_score_{args.version}.finite.txt").write_text("\n".join(finite) + "\n")

    if args.data == "trainval":
        scores = common.metrics_on(model, test, features)
        print(f"\nTEST (held out, n={len(test)})")
        for metric in common.METRICS:
            print(f"  {metric:>8}: {scores[metric]:.4f}")

    print(f"\n  -> {booster}  ({len(features)} features, {len(finite)} never missing)"
          f"\n  -> {MODEL_DIR / f'tauso_score_{args.version}.features.txt'}"
          f"\n  -> {MODEL_DIR / f'tauso_score_{args.version}.finite.txt'}")


if __name__ == "__main__":
    main()
