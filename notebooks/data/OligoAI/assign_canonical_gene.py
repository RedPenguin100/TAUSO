# notebooks/data/OligoAI/assign_canonical_gene.py
import pandas as pd
import warnings

from notebooks.consts import ORIGINAL_OLIGO_CSV, PROCESSED_OLIGO_CSV_GZ
from tauso.off_target.search import find_all_gene_off_targets
from tauso.data.consts import CANONICAL_GENE


def process_and_assign_genes(
    input_csv=ORIGINAL_OLIGO_CSV, output_csv=PROCESSED_OLIGO_CSV_GZ, genome="GRCh38"
):
    df = pd.read_csv(input_csv)

    missing_cols = {"rna_context", "target_gene"} - set(df.columns)
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    gene_to_context = df.dropna(subset=["target_gene", "rna_context"]).groupby("target_gene")["rna_context"].first()

    print(f"Resolving {len(gene_to_context)} unique target genes")

    gene_to_canonical = {}
    ambiguous_genes = {}

    for i, (target_gene, rna_seq) in enumerate(gene_to_context.items(), start=1):
        dna = rna_seq.upper().replace("U", "T")
        hits = find_all_gene_off_targets(dna, genome=genome, max_mismatches=0)
        genes = hits["gene_name"].dropna().unique().tolist() if not hits.empty else []

        if len(genes) == 0:
            gene_to_canonical[target_gene] = pd.NA
        elif len(genes) == 1:
            gene_to_canonical[target_gene] = genes[0]
        else:
            gene_to_canonical[target_gene] = ";".join(genes)
            ambiguous_genes[target_gene] = genes

        if i % 50 == 0 or i == len(gene_to_context):
            print(f"[{i}/{len(gene_to_context)}] genes processed")

    df[CANONICAL_GENE] = df["target_gene"].map(gene_to_canonical)
    df.to_csv(output_csv, index=False)
    print(f"Saved annotated dataframe → {output_csv}")

    if ambiguous_genes:
        warnings.warn(f"{len(ambiguous_genes)} target genes mapped to multiple canonical genes")
    else:
        print("No ambiguous gene mappings found 🎉")

    return df, ambiguous_genes


if __name__ == "__main__":
    # This block executes only when run from the command line (e.g., in CI)
    process_and_assign_genes()