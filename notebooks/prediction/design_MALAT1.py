"""Design ASOs against MALAT1 (nuclear lncRNA; looked up from the genome cache, no CDS so
target_region is `unknown`). Writes ranked / tox / feature CSVs to notebooks/prediction/output/vanilla/
(a separate set of CSVs per --chemistry preset). Use --first-n for a quick partial run."""

import argparse
import os
from pathlib import Path

from _chemistries import CHEMISTRIES, build_config

from tauso.aso_generation import default_config, design_asos, feature_details, summarize_design, tox_details

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chemistry", choices=list(CHEMISTRIES), default="moe_5_10_5", help="ASO chemistry preset")
    parser.add_argument("--first-n", type=int, default=None, help="featurize only the first N candidates (default: whole transcript)")
    parser.add_argument("--n-jobs", type=int, default=min(os.cpu_count() or 1, 24), help="worker processes")
    parser.add_argument("--out-dir", type=Path, default=HERE / "output" / "vanilla")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    config, aso_size, label = build_config(default_config, args.chemistry)

    ranked = design_asos(
        "MALAT1",
        aso_sizes=[aso_size],
        first_n=args.first_n,
        n_jobs=args.n_jobs,
        config=config,
    )

    designed = summarize_design(ranked)
    designed.to_csv(args.out_dir / f"MALAT1_{label}designed_asos.csv", index=False)
    tox_details(ranked).to_csv(args.out_dir / f"MALAT1_{label}tox.csv", index=False)
    feature_details(ranked).to_csv(args.out_dir / f"MALAT1_{label}features.csv", index=False)

    print(f"MALAT1 [{args.chemistry}]: {len(ranked)} candidate ASOs designed -> {args.out_dir}")
    print(designed.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
