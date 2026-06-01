"""Public per-gene hybridization features: score ASOs against a single target gene.

Thin entry points over the gene-chunk scoring orchestration — one for on-target (each
ASO vs its own canonical gene) and one for off-target against a fixed named gene. Each
returns the dataframe plus the names of the feature columns it added (one per cutoff).
"""

from ....data.consts import CANONICAL_GENE_NAME
from .gene_chunk_scoring import dispatch_gene_chunk_scores


def on_target_total_hybridization(aso_df, gene_to_data, cutoffs, n_jobs=1):
    """Score each oligo against its own canonical gene, one feature per cutoff."""
    return dispatch_gene_chunk_scores(
        aso_df=aso_df,
        gene_to_data=gene_to_data,
        target_genes=aso_df[CANONICAL_GENE_NAME].dropna().unique(),
        get_gene_fn=lambda row: row[CANONICAL_GENE_NAME],
        feature_name_fn=lambda c: f"on_target_total_hybridization_{c}",
        cutoffs=cutoffs,
        n_jobs=n_jobs,
    )


def off_target_single_gene_hybridization(aso_df, gene_name, gene_to_data, cutoffs, n_jobs=1):
    """Score each oligo against a single fixed gene, one feature per cutoff."""
    return dispatch_gene_chunk_scores(
        aso_df=aso_df,
        gene_to_data=gene_to_data,
        target_genes=[gene_name],
        get_gene_fn=lambda row: gene_name,
        feature_name_fn=lambda c: f"off_target_single_{gene_name}_c{c}",
        cutoffs=cutoffs,
        n_jobs=n_jobs,
    )
