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


def alignment_stats(input_csv=ORIGINAL_OLIGO_CSV_RAW, genome="GRCh38", seq_col="rna_sequence", threads=16):
    """Per patent target_gene: how often its ASOs align to that gene vs. the most popular hit."""
    df = pd.read_csv(input_csv)
    fasta = os.path.join(tempfile.gettempdir(), f"curate_{uuid.uuid4().hex}.fasta")
    with open(fasta, "w") as f:
        for seq in df[seq_col].dropna().unique():
            dna = seq.upper().replace("U", "T")
            f.write(f">{dna}\n{dna}\n")
    try:
        seq_to_genes = find_all_gene_off_targets_BULK(fasta, genome, threads)
    finally:
        os.path.exists(fasta) and os.remove(fasta)

    rows = []
    for gene, seqs in df.dropna(subset=["target_gene", seq_col]).groupby("target_gene")[seq_col].unique().items():
        hits = Counter()
        for seq in seqs:
            for g in seq_to_genes.get(seq.upper().replace("U", "T"), []):
                hits[g] += 1
        top = max(hits, key=hits.get) if hits else pd.NA
        rows.append({
            "original_target_gene": gene,
            "total_asos": len(seqs),
            "most_popular_alignment": top,
            "original_hit_count": hits.get(gene, 0),
            "popular_hit_count": hits.get(top, 0) if hits else 0,
        })
    return pd.DataFrame(rows)


def report(stats):
    disagree = stats[stats["original_target_gene"] != stats["most_popular_alignment"]].copy()
    disagree["original_below_half"] = disagree["original_hit_count"] < 0.5 * disagree["total_asos"]
    disagree["already_corrected"] = disagree["original_target_gene"].isin(MANUAL_CANONICAL_MAPPING)
    print(f"{len(disagree)} target_genes disagree with their most-popular alignment "
          f"({int(disagree['already_corrected'].sum())} already in MANUAL_CANONICAL_MAPPING):\n")
    print(disagree.sort_values("original_below_half", ascending=False).to_string(index=False))


if __name__ == "__main__":
    report(alignment_stats())
