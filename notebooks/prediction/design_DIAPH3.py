"""Generate features for DIAPH3 ASOs and rank them for OVCAR-8 / Lipofection / 100 nM.

design_asos tiles DIAPH3's pre-mRNA -- the genome lookup returns the 498 kb span with introns,
so intronic sites are tiled alongside exonic ones -- runs the tauso feature pipeline, scores with
the bundled model, and returns the candidates ranked best-first. Each candidate also gets
genome-wide off-target match counts (0/1/2 mismatches, the on-target DIAPH3 locus excluded) and a
liability block (off-target burden, rRNA binding, RNase-H1 fit, tox flags).

The chemistry is the library default: a 2'-MOE 5-10-5 gapmer on a full phosphorothioate backbone.
`--ps-pattern` swaps in a mixed PS/PO backbone, 19 characters of `*` (PS) and `d` (PO), one per
linkage, which feeds the mod_ps_* features.

  python notebooks/prediction/design_DIAPH3.py --first-n 2000
  python notebooks/prediction/design_DIAPH3.py                       # the full tiling
  python notebooks/prediction/design_DIAPH3.py --ps-pattern '***dddddddddddddd***'
"""
import argparse
import os
from pathlib import Path

import pandas as pd

from tauso.aso_generation import Transfection, default_config, design_asos, summarize_design, tox_details
from tauso.data.data import get_paths
from tauso.genome.read_human_genome import get_locus_to_data_dict
from tauso.off_target.search import count_offtarget_matches_bulk

HERE = Path(__file__).resolve().parent
OUT = HERE / "output" / "DIAPH3_ovcar8_lipo"

GENE = "DIAPH3"
GENOME = "GRCh38"
CELL_LINE = "OVCAR-8"
DELIVERY = Transfection.LIPOFECTION
DOSE_NM = 100

OFFTARGET_MAX_MM = 2  # count genome matches up to 2 mismatches (Bowtie -v caps at 3)
OFFTARGET_COLS = ["perfect_matches", "off_targets_1mm", "off_targets_2mm"]
LIABILITY_COLS = ["offtarget_transcriptome", "offtarget_genomewide", "offtarget_rrna", "liabilities"]


def build_config(ps_pattern):
    cfg = default_config()                          # 2'-MOE 5-10-5 gapmer, full-PS 20-mer
    cfg.cell_line = CELL_LINE
    cfg.transfection_method = DELIVERY
    cfg.volume = DOSE_NM                            # dose, nM
    cfg.organism_name = "human"
    cfg.standard_ps_pattern = ps_pattern
    return cfg


def add_offtarget_counts(summary):
    """Add per-ASO genome-wide Bowtie match counts at 0/1/2 mismatches, both strands, excluding
    hits inside the on-target DIAPH3 locus so the intended site and intragenic near-matches are not
    counted. Fills <NA> and prints a note if the Bowtie index is missing, rather than triggering a
    multi-GB index build (run `tauso setup-bowtie --genome GRCh38`)."""
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
    for col, i in zip(OFFTARGET_COLS, range(3)):
        summary[col] = summary["aso_sequence"].map(lambda s, i=i: counts[s][i])
    return summary


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--first-n", type=int, default=None, help="featurize only the first N candidates")
    ap.add_argument("--n-jobs", type=int, default=min(os.cpu_count() or 1, 24), help="worker processes")
    ap.add_argument("--ps-pattern", default=19 * "*",
                    help="backbone, one character per linkage: * is phosphorothioate, d phosphodiester")
    ap.add_argument("--tag", default=None, help="suffix for the output filenames")
    ap.add_argument("--no-offtargets", action="store_true", help="skip the genome-wide match counts")
    args = ap.parse_args()

    ranked = design_asos(GENE, config=build_config(args.ps_pattern), first_n=args.first_n, n_jobs=args.n_jobs)

    summary = summarize_design(ranked)
    if not args.no_offtargets:
        summary = add_offtarget_counts(summary)
    details = tox_details(ranked)                   # row-aligned to ranked -> summarize_design
    for col in LIABILITY_COLS:
        summary[col] = details[col].to_numpy()
    summary["chemistry"] = f"2'MOE 5-10-5 gapmer, backbone {args.ps_pattern}"
    summary["transfection_method"] = DELIVERY
    summary["dosage_nm"] = DOSE_NM
    summary["cell_line"] = CELL_LINE

    OUT.mkdir(parents=True, exist_ok=True)
    tag = f"_{args.tag}" if args.tag else ""
    ranked.to_parquet(OUT / f"DIAPH3_ovcar8_lipo{tag}_features.parquet", index=False)
    summary.to_csv(OUT / f"DIAPH3_ovcar8_lipo{tag}_ranked.csv", index=False)
    print(f"{len(ranked)} DIAPH3 ASOs featurized + ranked -> {OUT}")


if __name__ == "__main__":
    main()
