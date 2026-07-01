"""Design ASOs against the DIAPH3 plasmid construct (CDS 70..3652). Writes ranked / tox / feature
CSVs to notebooks/prediction/output/vanilla/. Use --first-n for a quick partial run."""

import argparse
import os
from pathlib import Path

from tauso.aso_generation import default_config, design_asos, feature_details, summarize_design, tox_details
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
    parser.add_argument("--out-dir", type=Path, default=HERE / "output" / "vanilla")
    parser.add_argument(
        "--ps-pattern",
        default=None,
        help="backbone pattern, one char per linkage ('*'=PS, 'd'=PO), length ASO_SIZE-1; default = full PS",
    )
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    config = default_config()
    if args.ps_pattern:
        if len(args.ps_pattern) != ASO_SIZE - 1:
            parser.error(f"--ps-pattern must be {ASO_SIZE - 1} chars for {ASO_SIZE}-mers, got {len(args.ps_pattern)}")
        config.standard_ps_pattern = args.ps_pattern

    _, seq = read_single_rna_fasta(str(FASTA))
    ranked = design_asos(
        "DIAPH3",
        gene_sequence=seq,
        cds_start=CDS_START,
        cds_end=CDS_END,
        aso_sizes=[ASO_SIZE],
        first_n=args.first_n,
        n_jobs=args.n_jobs,
        config=config,
    )

    designed = summarize_design(ranked)
    designed.to_csv(args.out_dir / "DIAPH3_designed_asos.csv", index=False)
    tox_details(ranked).to_csv(args.out_dir / "DIAPH3_tox.csv", index=False)
    feature_details(ranked).to_csv(args.out_dir / "DIAPH3_features.csv", index=False)

    print(f"DIAPH3: {len(ranked)} candidate ASOs designed -> {args.out_dir}")
    print(designed.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
