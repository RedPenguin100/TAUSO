"""On-target multiplicity: log effective number of on-target sites per ASO.

A thin entry point over the per-gene scan (each ASO vs its own canonical gene), deriving the log
effective site count from the site-resolved stats (see
off_target_specific_gene.log_number_of_sites_from_stats). The pipeline emits this column via
``add_on_target_site_features``; this entry point runs its own scan and backs the regression test.
"""

from ....data.consts import CANONICAL_GENE_NAME
from .gene_chunk_scoring import emit_site_columns, scan_gene_sites
from .off_target_specific_gene import log_number_of_sites_from_stats


def on_target_log_number_of_sites(aso_df, gene_to_data, cutoffs, n_jobs=1):
    """Log effective number of on-target sites per ASO, one column per cutoff.

    Each oligo is scored against its own canonical gene. Returns (aso_df, [feature_names]) with
    columns ``on_target_log_number_of_sites_{cutoff}``. ASOs with no gene or no qualifying hit get 0.
    """
    scan = scan_gene_sites(
        aso_df=aso_df,
        gene_to_data=gene_to_data,
        target_genes=aso_df[CANONICAL_GENE_NAME].dropna().unique(),
        get_gene_fn=lambda row: row[CANONICAL_GENE_NAME],
        cutoffs=cutoffs,
        n_jobs=n_jobs,
    )
    return emit_site_columns(
        aso_df,
        scan,
        cutoffs,
        derivations=[(lambda c: f"on_target_log_number_of_sites_{c}", log_number_of_sites_from_stats)],
    )
