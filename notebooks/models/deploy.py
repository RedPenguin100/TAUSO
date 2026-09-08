"""Deploy -- train and save the final shipped model.

Trains the config in deploy_parameters.json on train+val and writes the booster, the feature
list the package scores with, and the list of those features that are never missing.

One of --use-calculated / --use-downloaded says where the features come from. They are
different feature sets whenever the pipeline has moved on since the cache was published, so
the choice is the caller's rather than whichever happens to be on disk.

  python notebooks/models/deploy.py --use-calculated          # what `calculate_features` wrote
  python notebooks/models/deploy.py --use-downloaded          # the published cache, fetched if absent
  python notebooks/models/deploy.py --use-calculated --med    # the MED search's parameters
  python notebooks/models/deploy.py --use-calculated --data all

The booster is ~100 MB and stays out of git; copy it to <data_dir>/models/ to score with it.
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # repo root, for the notebooks.* imports
from notebooks.models import common

from tauso.inference.predict import DEFAULT_VERSION, MODEL_DIR, MODEL_FILES

CONFIGS = json.loads((Path(__file__).parent / "deploy_parameters.json").read_text())


def shard_dir():
    """Where the feature pipeline writes its per-feature shards."""
    from notebooks.features.feature_extraction import _get_saved_features_dir

    from tauso.populate.feature_cache import loose_shard_dir

    return Path(loose_shard_dir(_get_saved_features_dir("oligo")))


def current_pipeline_features():
    """Names the feature pipeline writes a shard for, i.e. what it still computes."""
    directory = shard_dir()
    if not directory.is_dir():
        return set()
    return {p.stem for p in directory.iterdir() if p.suffix in (".parquet", ".csv")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--data",
        choices=["trainval", "all"],
        default="trainval",
        help="train on train+val (default, the shipped model) or all data",
    )
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--version", default=DEFAULT_VERSION, help="model version to write")
    config = ap.add_mutually_exclusive_group()
    config.add_argument("--low", action="store_true", help="the LOW search's parameters (default)")
    config.add_argument("--med", action="store_true", help="the MED search's parameters")
    source = ap.add_mutually_exclusive_group(required=True)
    source.add_argument("--use-calculated", action="store_true", help="train on the features this machine computed")
    source.add_argument(
        "--use-downloaded", action="store_true", help="train on the published cache, fetching it if absent"
    )
    args = ap.parse_args()

    config_name = "med" if args.med else "low"
    spec = CONFIGS[config_name]

    if args.use_calculated and not current_pipeline_features():
        sys.exit(
            f"No calculated features in {shard_dir()}.\n"
            "Compute them first:\n"
            "  python -m notebooks.features.calculate_features --dataset oligo --cpus $(nproc)\n"
            "or train on the published set with --use-downloaded."
        )

    df, features = common.load_dataset(use_cache=args.use_downloaded)
    print(f"features from the {'published cache' if args.use_downloaded else 'feature pipeline'}", flush=True)

    if args.use_calculated:
        # Nothing may come from the cache: it still holds columns the pipeline has stopped
        # producing, and a model trained on one of those cannot be scored by a feature run.
        produced = current_pipeline_features()
        strays = sorted(set(features) - produced)
        if strays:
            sys.exit(f"{len(strays)} features are not among the calculated shards: {strays}")

    # A column with one value carries no split. All-NaN is the degenerate case, and it means
    # the data behind the feature was never built rather than the feature being uninformative.
    constant = [f for f in features if df[f].nunique(dropna=True) <= 1]
    empty = [f for f in constant if df[f].isna().all()]
    if empty:
        print(f"\n!! dropping {len(empty)} features that are NaN for every row: {empty}")
        print("!! the data they are built from is missing; they will be absent from the model\n", flush=True)
    if len(constant) > len(empty):
        print(
            f"dropping {len(constant) - len(empty)} constant features: {[f for f in constant if f not in set(empty)]}",
            flush=True,
        )
    features = [f for f in features if f not in set(constant)]
    trv, test = common.split(df)
    train_df = trv if args.data == "trainval" else df

    params, rounds, variant = spec["params"], spec["num_boost_round"], spec["variant"]
    print(
        f"deploy {args.version} ({variant}, {config_name}): train on {args.data} ({len(train_df)} rows), "
        f"seed {args.seed}, {rounds} rounds, {len(features)} feats",
        flush=True,
    )

    model = common.train(train_df, features, variant, params, rounds, seed=args.seed)

    booster = common.RESULTS_DIR / f"{Path(MODEL_FILES[args.version]['filename']).stem}_{config_name}.json"
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

    print(
        f"\n  -> {booster}  ({len(features)} features, {len(finite)} never missing)"
        f"\n  -> {MODEL_DIR / f'tauso_score_{args.version}.features.txt'}"
        f"\n  -> {MODEL_DIR / f'tauso_score_{args.version}.finite.txt'}"
    )


if __name__ == "__main__":
    main()
