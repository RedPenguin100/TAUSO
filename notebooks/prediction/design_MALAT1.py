"""Generate v15 features for MALAT1 ASOs and rank them for A549 / Lipofection / 100 nM.

MALAT1 is a nuclear lncRNA; its sequence is looked up from the genome cache (no CDS, so the
target region is `unknown`). design_asos tiles it into 20-mers, runs the tauso feature pipeline
(calculate_features), scores with the bundled model, and returns the candidates ranked best-first.
Each candidate also gets genome-wide off-target match counts (0/1/2 mismatches, on-target MALAT1
locus excluded) and a liability block (off-target burden, rRNA binding, RNase-H1 fit, tox flags).

  python notebooks/prediction/design_MALAT1.py
  python notebooks/prediction/design_MALAT1.py --first-n 200   # quick partial run
"""
import argparse
import os
from pathlib import Path

import pandas as pd

from tauso.aso_generation import Transfection, default_config, design_asos, design_details, summarize_design
from tauso.data.data import get_paths
from tauso.genome.read_human_genome import get_locus_to_data_dict
from tauso.off_target.search import count_offtarget_matches_bulk

HERE = Path(__file__).resolve().parent
OUT = HERE / "output" / "MALAT1_A549_lipo"

GENE = "MALAT1"
GENOME = "GRCh38"
CELL_LINE = "A549"
DELIVERY = Transfection.LIPOFECTION
DOSE_NM = 100
CHEMISTRY = "2'MOE 5-10-5 gapmer, full PS"

OFFTARGET_MAX_MM = 2  # count genome matches up to 2 mismatches (Bowtie -v caps at 3)
OFFTARGET_COLS = ["perfect_matches", "off_targets_1mm", "off_targets_2mm"]
LIABILITY_COLS = ["offtarget_transcriptome", "offtarget_genomewide", "offtarget_rrna",
                  "rnaseh1_cleavage_fit", "liabilities"]


def build_config():
    cfg = default_config()                         # 2'-MOE 5-10-5 gapmer, full-PS 20-mer
    cfg.cell_line = CELL_LINE
    cfg.transfection_method = DELIVERY
    cfg.volume = DOSE_NM                            # dose, nM
    cfg.organism_name = "human"
    return cfg


def add_offtarget_counts(summary):
    """Add per-ASO off-target counts: `perfect_matches` (0 mismatches), `off_targets_1mm`, and
    `off_targets_2mm` -- genome-wide Bowtie match counts at 0/1/2 mismatches (both strands, one pass),
    EXCLUDING any hit inside the on-target MALAT1 locus so the intended site and intragenic near-matches
    are not counted. If the Bowtie index is missing the columns are filled with <NA> and a note printed,
    rather than triggering a multi-GB index build (run `tauso setup-bowtie --genome GRCh38`)."""
    sentinel = Path(get_paths(GENOME)["fasta"]).parent / f"{GENOME}_bowtie_index" / "SUCCESS"
    if not sentinel.exists():
        print(f"Bowtie index for {GENOME} not found; skipping off-target counts "
              f"(run `tauso setup-bowtie --genome {GENOME}`).")
        for col in OFFTARGET_COLS:
            summary[col] = pd.NA
        return summary

    g = get_locus_to_data_dict(include_introns=True, gene_subset=[GENE], genome=GENOME)[GENE]
    counts = count_offtarget_matches_bulk(
        summary["aso_sequence"].tolist(), genome=GENOME,
        max_mismatches=OFFTARGET_MAX_MM, exclude_regions=[(g.chrom, g.gene_start, g.gene_end)],
    )
    summary["perfect_matches"] = summary["aso_sequence"].map(lambda s: counts[s][0])
    summary["off_targets_1mm"] = summary["aso_sequence"].map(lambda s: counts[s][1])
    summary["off_targets_2mm"] = summary["aso_sequence"].map(lambda s: counts[s][2])
    return summary


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--first-n", type=int, default=None, help="featurize only the first N candidates")
    ap.add_argument("--n-jobs", type=int, default=min(os.cpu_count() or 1, 24), help="worker processes")
    args = ap.parse_args()

    ranked = design_asos(GENE, config=build_config(), first_n=args.first_n, n_jobs=args.n_jobs)

    summary = summarize_design(ranked)
    summary = add_offtarget_counts(summary)
    details = design_details(ranked)                       # row-aligned to ranked -> summarize_design
    for col in LIABILITY_COLS:
        summary[col] = details[col].to_numpy()
    summary["chemistry"] = CHEMISTRY
    summary["transfection_method"] = DELIVERY
    summary["dosage_nm"] = DOSE_NM

    OUT.mkdir(parents=True, exist_ok=True)
    ranked.to_parquet(OUT / "MALAT1_A549_lipo_features.parquet", index=False)
    summary.to_csv(OUT / "MALAT1_A549_lipo_ranked.csv", index=False)
    print(f"{len(ranked)} MALAT1 ASOs featurized + ranked (off-target counts + liabilities) -> {OUT}")


if __name__ == "__main__":
    main()
