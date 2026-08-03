"""Diagnostic: align every ASO to the genome and flag patents whose target_gene disagrees with the
gene the ASOs actually align to. This is how MANUAL_CANONICAL_MAPPING (gene_corrections.py) was
derived; it is NOT part of the per-run pipeline (2_assign_canonical_gene applies the frozen map).

Run: python notebooks/data/OligoAI/curate_gene_labels.py   (needs the bowtie index; tauso setup-all)
"""
import ast
import os
import sys
import tempfile
import uuid
from collections import Counter
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from notebooks.consts import ORIGINAL_OLIGO_CSV_RAW
from notebooks.data.OligoAI.gene_corrections import MANUAL_CANONICAL_MAPPING
from tauso.off_target.search import find_all_gene_off_targets_BULK


def alignment_stats(input_csv=ORIGINAL_OLIGO_CSV_RAW, genome="GRCh38", seq_col="aso_sequence_5_to_3", threads=16):
    """Per patent label: how often its ASOs align to that gene vs. the most popular hit.

    Counts come in both units: distinct sequences (`_asos`) and measurements (`_measurements`).
    """
    df = pd.read_csv(input_csv).dropna(subset=[seq_col])
    # Patents that leave target_gene blank still name the target in target_mrna, as 2_assign_canonical_gene
    # assumes. Falling back unconditionally is what lets this diagnostic find cohorts not yet in the map.
    fallback = df["target_mrna"] if "target_mrna" in df else pd.Series(pd.NA, index=df.index)
    label = df["target_gene"].fillna(fallback)
    source = df["target_gene"].notna().map({True: "target_gene", False: "target_mrna"})
    df, label, source = df[label.notna()], label.dropna(), source[label.notna()]
    dna = df[seq_col].str.upper().str.replace("U", "T", regex=False)

    fasta = os.path.join(tempfile.gettempdir(), f"curate_{uuid.uuid4().hex}.fasta")
    with open(fasta, "w") as f:
        for seq in dna.unique():
            f.write(f">{seq}\n{seq}\n")
    try:
        seq_to_genes = find_all_gene_off_targets_BULK(fasta, genome, threads)
    finally:
        os.path.exists(fasta) and os.remove(fasta)

    rows = []
    for gene, group in dna.groupby(label):
        asos, measurements = Counter(), Counter()
        for seq, n in group.value_counts().items():
            for g in seq_to_genes.get(seq, []):
                asos[g] += 1
                measurements[g] += n
        top = max(asos, key=asos.get) if asos else pd.NA
        rows.append({
            "original_target_gene": gene,
            "label_from": source[group.index].mode().iat[0],
            "total_asos": group.nunique(),
            "total_measurements": len(group),
            "most_popular_alignment": top,
            "original_hit_asos": asos[gene],
            "original_hit_measurements": measurements[gene],
            "popular_hit_asos": asos[top] if asos else 0,
            "popular_hit_measurements": measurements[top] if asos else 0,
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
