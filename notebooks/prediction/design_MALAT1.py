"""Generate v15 features for MALAT1 ASOs and rank them for A549 / Lipofection / 100 nM.

MALAT1 is a nuclear lncRNA; its sequence is looked up from the genome cache (no CDS, so the
target region is `unknown`). design_asos tiles it into 20-mers, runs the tauso feature pipeline
(calculate_features), scores with the bundled model, and returns the candidates ranked best-first.

  python notebooks/prediction/design_MALAT1.py
  python notebooks/prediction/design_MALAT1.py --first-n 200   # quick partial run
"""
import argparse
import os
from pathlib import Path

from tauso.aso_generation import Transfection, default_config, design_asos, summarize_design

HERE = Path(__file__).resolve().parent
OUT = HERE / "output" / "MALAT1_A549_lipo"


def build_config():
    cfg = default_config()                         # 2'-MOE 5-10-5 gapmer, full-PS 20-mer
    cfg.cell_line = "A549"
    cfg.transfection_method = Transfection.LIPOFECTION
    cfg.volume = 100                               # dose, nM
    cfg.organism_name = "human"
    return cfg


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--first-n", type=int, default=None, help="featurize only the first N candidates")
    ap.add_argument("--n-jobs", type=int, default=min(os.cpu_count() or 1, 24), help="worker processes")
    args = ap.parse_args()

    ranked = design_asos("MALAT1", config=build_config(), first_n=args.first_n, n_jobs=args.n_jobs)

    OUT.mkdir(parents=True, exist_ok=True)
    ranked.to_parquet(OUT / "MALAT1_A549_lipo_features.parquet", index=False)
    summarize_design(ranked).to_csv(OUT / "MALAT1_A549_lipo_ranked.csv", index=False)
    print(f"{len(ranked)} MALAT1 ASOs featurized + ranked -> {OUT}")


if __name__ == "__main__":
    main()
