"""Per-gene hybridization features, all derived from a RIsearch scan.

Each entry point runs the per-gene scan (gene_chunk_scoring.scan_gene_sites) and derives its
columns from the site-resolved stats:

- ``on_target_total_hybridization``  -- each ASO vs its own canonical gene, Sum exp(-energy/RT)
- ``on_target_log_number_of_sites``  -- same scan, log effective number of sites (multiplicity)
- ``off_target_single_gene_hybridization`` -- each ASO vs a fixed named gene, Sum exp(-energy/RT)

``add_on_target_site_features`` emits BOTH on-target features together; it is what the
feature pipeline calls. The individual on-target entry points each run their own scan and back the
regression tests.
"""

import numpy as np

from ....data.consts import CANONICAL_GENE_NAME
from ..fast_hybridization import RT_KCAL_MOL
from .gene_chunk_scoring import SiteStats, emit_site_columns, scan_gene_sites


def total_hybridization_from_stats(stats: SiteStats) -> float:
    """Total hybridization = the Boltzmann occupancy sum over the ASO's sites."""
    return stats.sum_exp


def log_number_of_sites_from_stats(stats: SiteStats):
    """Log effective number of sites (target multiplicity): 0 for a single dominant site, growing
    when several comparable sites share the binding. None when the ASO has no occupancy."""
    if stats.sum_exp <= 0:
        return None
    log_eff = np.log(max(stats.sum_exp, 1e-300)) - (-stats.min_energy / RT_KCAL_MOL)
    return log_eff if np.isfinite(log_eff) else None


def _scan_own_gene(aso_df, gene_to_data, cutoffs, n_jobs):
    """Scan each ASO against its own canonical gene."""
    return scan_gene_sites(
        aso_df=aso_df,
        gene_to_data=gene_to_data,
        target_genes=aso_df[CANONICAL_GENE_NAME].dropna().unique(),
        get_gene_fn=lambda row: row[CANONICAL_GENE_NAME],
        cutoffs=cutoffs,
        n_jobs=n_jobs,
    )


def add_on_target_site_features(aso_df, gene_to_data, cutoffs, n_jobs=1):
    """Scores each ASO vs its own canonical gene -> BOTH on-target features per
    cutoff: ``on_target_total_hybridization_{c}`` and ``on_target_log_number_of_sites_{c}``.

    Returns (aso_df, [feature_names]).
    """
    scan = _scan_own_gene(aso_df, gene_to_data, cutoffs, n_jobs)
    return emit_site_columns(
        aso_df,
        scan,
        cutoffs,
        derivations=[
            (lambda c: f"on_target_total_hybridization_{c}", total_hybridization_from_stats),
            (lambda c: f"on_target_log_number_of_sites_{c}", log_number_of_sites_from_stats),
        ],
    )


def on_target_total_hybridization(aso_df, gene_to_data, cutoffs, n_jobs=1):
    """Score each oligo against its own canonical gene: total hybridization, one feature per cutoff."""
    scan = _scan_own_gene(aso_df, gene_to_data, cutoffs, n_jobs)
    return emit_site_columns(
        aso_df,
        scan,
        cutoffs,
        derivations=[(lambda c: f"on_target_total_hybridization_{c}", total_hybridization_from_stats)],
    )


def off_target_single_gene_hybridization(aso_df, gene_name, gene_to_data, cutoffs, n_jobs=1):
    """Score each oligo against a single fixed gene: total hybridization, one feature per cutoff."""
    scan = scan_gene_sites(
        aso_df=aso_df,
        gene_to_data=gene_to_data,
        target_genes=[gene_name],
        get_gene_fn=lambda row: gene_name,
        cutoffs=cutoffs,
        n_jobs=n_jobs,
    )
    return emit_site_columns(
        aso_df,
        scan,
        cutoffs,
        derivations=[(lambda c: f"off_target_single_{gene_name}_c{c}", total_hybridization_from_stats)],
    )
