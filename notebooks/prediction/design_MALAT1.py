"""Design ASOs against MALAT1 (nuclear lncRNA).

MALAT1 is a single-exon, nuclear-retained long non-coding RNA -- a classic gapmer target -- so
it has no coding region. The target transcript is looked up from the genome cache (no FASTA and
no CDS annotation); consequently every candidate's ``target_region`` is reported as ``unknown``
(there is no 5'UTR / CDS / 3'UTR to assign).

Tiles every 20-mer across the transcript, scores each with the bundled ``tauso_score_v1`` model,
and writes two CSVs to ``notebooks/prediction/output/`` (same schema as design_DIAPH3.py):

  * ``MALAT1_designed_asos.csv``  -- ranked sequences.
  * ``MALAT1_safety_detail.csv``  -- per-ASO off-target / rRNA / RNase H1 fit and toxicity flags.

Run:  python notebooks/prediction/design_MALAT1.py            # full transcript (~8.8 kb)
      python notebooks/prediction/design_MALAT1.py --first-n 20   # quick smoke run
"""

import argparse
import os
from pathlib import Path

from tauso.aso_generation import design_asos, design_details, summarize_design

HERE = Path(__file__).resolve().parent
ASO_SIZE = 20


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first-n", type=int, default=None, help="featurize only the first N candidates (default: whole transcript)")
    parser.add_argument("--n-jobs", type=int, default=min(os.cpu_count() or 1, 24), help="worker processes")
    parser.add_argument("--out-dir", type=Path, default=HERE / "output")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    ranked = design_asos(
        "MALAT1",
        aso_sizes=[ASO_SIZE],
        first_n=args.first_n,
        n_jobs=args.n_jobs,
    )

    designed = summarize_design(ranked)
    safety = design_details(ranked)
    designed.to_csv(args.out_dir / "MALAT1_designed_asos.csv", index=False)
    safety.to_csv(args.out_dir / "MALAT1_safety_detail.csv", index=False)

    print(f"MALAT1: {len(ranked)} candidate ASOs designed -> {args.out_dir}")
    print(designed.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
