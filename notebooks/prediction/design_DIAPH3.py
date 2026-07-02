"""Design ASOs against the DIAPH3 plasmid construct (CDS 70..3652). Writes ranked / tox / feature
CSVs to notebooks/prediction/output/vanilla/ (a separate set of CSVs per --chemistry preset). Use
--first-n for a quick partial run."""

import argparse
import os
from pathlib import Path

from _chemistries import CHEMISTRIES, build_config

from tauso.aso_generation import default_config, design_asos, feature_details, summarize_design, tox_details
from tauso.genome.fasta import read_single_rna_fasta

HERE = Path(__file__).resolve().parent
FASTA = HERE / "DIAPH3_plasmid_sequence.fasta"

# DIAPH3 plasmid ORF / CDS = [70, 3652) (0-based, half-open) -> 5'UTR 70 / CDS 3582 / 3'UTR 133 nt.
CDS_START, CDS_END = 70, 3652


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chemistry", choices=list(CHEMISTRIES), default="moe_5_10_5", help="ASO chemistry preset")
    parser.add_argument("--first-n", type=int, default=None, help="featurize only the first N candidates (default: whole transcript)")
    parser.add_argument("--n-jobs", type=int, default=min(os.cpu_count() or 1, 24), help="worker processes")
    parser.add_argument("--out-dir", type=Path, default=HERE / "output")
    args = parser.parse_args()

    config, aso_size, subdir = build_config(default_config, args.chemistry)
    out_dir = args.out_dir / subdir
    out_dir.mkdir(parents=True, exist_ok=True)

    _, seq = read_single_rna_fasta(str(FASTA))
    ranked = design_asos(
        "DIAPH3",
        gene_sequence=seq,
        cds_start=CDS_START,
        cds_end=CDS_END,
        aso_sizes=[aso_size],
        first_n=args.first_n,
        n_jobs=args.n_jobs,
        config=config,
    )

    designed = summarize_design(ranked)
    designed.to_csv(out_dir / "DIAPH3_designed_asos.csv", index=False)
    tox_details(ranked).to_csv(out_dir / "DIAPH3_tox.csv", index=False)
    feature_details(ranked).to_csv(out_dir / "DIAPH3_features.csv", index=False)

    print(f"DIAPH3 [{args.chemistry}]: {len(ranked)} candidate ASOs designed -> {out_dir}")
    print(designed.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
