"""Diagnostic: align every ASO to the genome and flag patents whose target_gene disagrees with the
gene the ASOs actually align to. This is how MANUAL_CANONICAL_MAPPING (gene_corrections.py) was
derived; it is NOT part of the per-run pipeline (2_assign_canonical_gene applies the frozen map).

Run: python notebooks/data/OligoAI/curate_gene_labels.py   (needs the bowtie index; tauso setup-all)
"""
import ast
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from notebooks.consts import ORIGINAL_OLIGO_CSV_RAW
from notebooks.data.OligoAI.gene_corrections import MANUAL_CANONICAL_MAPPING
from tauso.off_target.search import find_all_gene_off_targets_bulk_sequences


def alignment_stats(input_csv=ORIGINAL_OLIGO_CSV_RAW, genome="GRCh38", threads=16):
    """Per patent label: how often its ASOs align to that gene vs. the most popular hit.

    Counts come in both units: distinct sequences (`_asos`) and measurements (`_measurements`).
    """
    SEQ_COL = "aso_sequence_5_to_3"

    df = pd.read_csv(input_csv)
    missing = [c for c in (SEQ_COL, "target_gene", "target_mrna") if c not in df.columns]
    if missing:
        raise KeyError(f"{input_csv} is missing required column(s): {missing}")
    df = df.dropna(subset=[SEQ_COL]).copy()

    # Patents that leave target_gene blank still name the target in target_mrna, as 2_assign_canonical_gene
    # assumes. Falling back unconditionally is what lets this diagnostic find cohorts not yet in the map.
    df["label"] = df["target_gene"].fillna(df["target_mrna"])
    df["label_from"] = np.where(df["target_gene"].notna(), "target_gene", "target_mrna")
    df = df.dropna(subset=["label"])
    df["aso_as_dna"] = df[SEQ_COL].str.upper().str.replace("U", "T", regex=False)

    seq_to_genes = find_all_gene_off_targets_bulk_sequences(df["aso_as_dna"].unique(), genome, threads)

    rows = []
    for gene, group in df.groupby("label"):
        asos, measurements = Counter(), Counter()
        for seq, n in group["aso_as_dna"].value_counts().items():
            for g in seq_to_genes.get(seq, []):
                asos[g] += 1
                measurements[g] += n
        top = max(asos, key=asos.get) if asos else pd.NA
        rows.append({
            "original_target_gene": gene,
            "label_from": group["label_from"].iat[0],
            "total_asos": group["aso_as_dna"].nunique(),
            "total_measurements": len(group),
            "most_popular_alignment": top,
            "original_hit_asos": asos[gene],
            "original_hit_measurements": measurements[gene],
            "popular_hit_asos": asos[top],
            "popular_hit_measurements": measurements[top],
        })
    return pd.DataFrame(rows)


def report(stats):
    disagree = stats[stats["original_target_gene"] != stats["most_popular_alignment"]].copy()
    disagree["original_below_half"] = disagree["original_hit_asos"] < 0.5 * disagree["total_asos"]
    disagree["already_corrected"] = disagree["original_target_gene"].isin(MANUAL_CANONICAL_MAPPING)
    print(f"{len(disagree)} target_genes disagree with their most-popular alignment "
          f"({int(disagree['already_corrected'].sum())} already in MANUAL_CANONICAL_MAPPING):\n")
    print(disagree.sort_values("original_below_half", ascending=False).to_string(index=False))


if __name__ == "__main__":
    report(alignment_stats())
