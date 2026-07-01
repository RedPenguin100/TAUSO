"""Design ASOs against the DIAPH3 plasmid construct.

Tiles every 20-mer across the DIAPH3 plasmid sequence, scores each with the bundled
`tauso_score_v1` efficacy model, and writes two CSVs to ``notebooks/prediction/output/``:

  * ``DIAPH3_designed_asos.csv``  -- ranked sequences (rank, gene, aso_sequence, length,
    target_start, target_region, tauso_score_v1). Higher score = better predicted knockdown
    (a within-experiment rank, not a percent-inhibition value).
  * ``DIAPH3_safety_detail.csv``  -- per-ASO off-target / rRNA / RNase H1 fit and toxicity
    liabilities (CpG/immune, G-quadruplex / G-run), with flags and an implications note.

Run:  python notebooks/prediction/design_DIAPH3.py            # full transcript
      python notebooks/prediction/design_DIAPH3.py --first-n 20   # quick smoke run
"""

import argparse
import os
from pathlib import Path

from tauso.aso_generation import design_asos, design_details, summarize_design
from tauso.genome.fasta import read_single_rna_fasta

HERE = Path(__file__).resolve().parent
FASTA = HERE / "DIAPH3_plasmid_sequence.fasta"

# DIAPH3 plasmid ORF / CDS = [70, 3652) (0-based, half-open) -> 5'UTR 70 / CDS 3582 / 3'UTR 133 nt.
CDS_START, CDS_END = 70, 3652
ASO_SIZE = 20


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first-n", type=int, default=None, help="featurize only the first N candidates (default: whole transcript)")
    parser.add_argument("--n-jobs", type=int, default=min(os.cpu_count() or 1, 24), help="worker processes")
    parser.add_argument("--out-dir", type=Path, default=HERE / "output")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    _, seq = read_single_rna_fasta(str(FASTA))
    ranked = design_asos(
        "DIAPH3",
        gene_sequence=seq,
        cds_start=CDS_START,
        cds_end=CDS_END,
        aso_sizes=[ASO_SIZE],
        first_n=args.first_n,
        n_jobs=args.n_jobs,
    )

    designed = summarize_design(ranked)
    safety = design_details(ranked)
    designed.to_csv(args.out_dir / "DIAPH3_designed_asos.csv", index=False)
    safety.to_csv(args.out_dir / "DIAPH3_safety_detail.csv", index=False)

    print(f"DIAPH3: {len(ranked)} candidate ASOs designed -> {args.out_dir}")
    print(designed.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
