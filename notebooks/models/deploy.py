"""Deploy -- train and save the final shipped model.

Trains the best clean_exp champion (top cross-validated exp_med on the v15 feature set) on
train+val and writes the deployable booster + its feature list. This is the model we ship.

  python notebooks/models/deploy.py                 # train on train+val, seed 1  (the final model)
  python notebooks/models/deploy.py --data all      # train on all data (train+val+test)
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # repo root, for the notebooks.* imports
from notebooks.models import champions, common

CHAMPION = "clean_exp_exp_med"     # the best clean_exp config
LABEL = "best_clean_exp"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", choices=["trainval", "all"], default="trainval",
                    help="train on train+val (default, the shipped model) or all data")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    df, features = common.load_dataset()
    features = [f for f in features if df[f].nunique(dropna=True) > 1]   # drop zero-variance -> the deployed 485-feature set
    trv, _ = common.split(df)
    train_df = trv if args.data == "trainval" else df

    params, rounds, variant = champions.config(CHAMPION)
    print(f"deploy 'best clean_exp' ({variant}): train on {args.data} "
          f"({len(train_df)} rows), seed {args.seed}, {rounds} rounds, {len(features)} feats", flush=True)

    model = common.train(train_df, features, variant, params, rounds, seed=args.seed)

    stem = common.RESULTS_DIR / f"deploy_model_{LABEL}_{args.data}_seed{args.seed}"
    model.save_model(str(stem.with_suffix(".json")))
    stem.with_suffix(".features.txt").write_text("\n".join(features) + "\n")
    print(f"deployed 'best clean_exp' ({variant}) on {args.data} "
          f"({len(train_df)} rows), seed {args.seed}, {len(features)} feats\n"
          f"  -> {stem.with_suffix('.json')}\n"
          f"  -> {stem.with_suffix('.features.txt')}")


if __name__ == "__main__":
    main()
