from tauso.off_target.search import find_all_gene_off_targets
import pandas as pd
import warnings

# =========================
# Config
# =========================

INPUT_CSV = "aso_inhibitions_21_08_25_incl_context_w_flank_100_df.csv.gz"
OUTPUT_CSV = "aso_inhibitions_with_canonical_gene.csv"
AMBIGUOUS_CSV = "ambiguous_target_gene_mappings.csv"

CONTEXT_COL = "rna_context"
TARGET_GENE_COL = "target_gene"
CANONICAL_COL = "Canonical Gene Name"

GENOME = "GRCh38"
MAX_MISMATCHES = 0   # for >200bp context, exact match only


# =========================
# Load data
# =========================

df = pd.read_csv(INPUT_CSV)

missing_cols = {CONTEXT_COL, TARGET_GENE_COL} - set(df.columns)
if missing_cols:
    raise ValueError(f"Missing required columns: {missing_cols}")


# =========================
# One context per target gene
# =========================

gene_to_context = (
    df
    .dropna(subset=[TARGET_GENE_COL, CONTEXT_COL])
    .groupby(TARGET_GENE_COL)[CONTEXT_COL]
    .first()
)

print(f"Resolving {len(gene_to_context)} unique target genes")


# =========================
# Resolve canonical genes
# =========================

gene_to_canonical = {}
ambiguous_genes = {}

for i, (target_gene, rna_seq) in enumerate(gene_to_context.items(), start=1):
    dna = rna_seq.upper().replace("U", "T")

    hits = find_all_gene_off_targets(
        dna,
        genome=GENOME,
        max_mismatches=MAX_MISMATCHES
    )

    genes = (
        hits["gene_name"]
        .dropna()
        .unique()
        .tolist()
        if not hits.empty
        else []
    )

    if len(genes) == 0:
        gene_to_canonical[target_gene] = pd.NA

    elif len(genes) == 1:
        gene_to_canonical[target_gene] = genes[0]

    else:
        gene_to_canonical[target_gene] = ";".join(genes)
        ambiguous_genes[target_gene] = genes

    if i % 50 == 0 or i == len(gene_to_context):
        print(f"[{i}/{len(gene_to_context)}] genes processed")


# =========================
# Map back to full dataframe
# =========================

df[CANONICAL_COL] = df[TARGET_GENE_COL].map(gene_to_canonical)


# =========================
# Save results
# =========================

df.to_csv(OUTPUT_CSV, index=False)
print(f"Saved annotated dataframe → {OUTPUT_CSV}")


# =========================
# Save ambiguous mappings
# =========================

if ambiguous_genes:
    warnings.warn(
        f"{len(ambiguous_genes)} target genes mapped to multiple canonical genes"
    )

    amb_df = pd.DataFrame(
        [
            {
                TARGET_GENE_COL: tg,
                "matched_canonical_genes": ";".join(genes)
            }
            for tg, genes in ambiguous_genes.items()
        ]
    )

    amb_df.to_csv(AMBIGUOUS_CSV, index=False)
    print(f"Saved ambiguous mappings → {AMBIGUOUS_CSV}")
else:
    print("No ambiguous gene mappings found 🎉")