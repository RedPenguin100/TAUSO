"""Per-gene hybridization features, all derived from a RIsearch scan (gene_chunk_scoring).

- ``on_target_total_hybridization``  -- each ASO vs its own canonical gene: Σ exp(-energy/RT)
- ``on_target_log_number_of_sites``  -- same scan: log effective number of sites (multiplicity)
- ``off_target_single_gene_hybridization`` -- each ASO vs a fixed named gene: Σ exp(-energy/RT)

``add_on_target_site_features`` emits both on-target features together (the pipeline entry point);
the individual entry points run their own scan and back the regression tests.
"""

import numpy as np

from ....data.consts import CANONICAL_GENE_NAME
from ....pandas_utils import add_columns
from . import RT_KCAL_MOL
from .gene_chunk_scoring import scan_gene_sites


def total_hybridization(sites, cutoff):
    """The Boltzmann occupancy summed over an ASO's sites. No hit is no occupancy."""
    return sites[cutoff, "sum_exp"].fillna(0.0)


def log_number_of_sites(sites, cutoff):
    """Log effective number of sites: 0 for a single dominant site, growing when
    several comparable sites share the binding.

    An ASO with no occupancy has no multiplicity to speak of and is given 0.0,
    which is also what one dominant site scores.
    """
    sum_exp = sites[cutoff, "sum_exp"]
    log_eff = np.log(sum_exp.clip(lower=1e-300)) - (-sites[cutoff, "min_energy"] / RT_KCAL_MOL)
    return log_eff.where(np.isfinite(log_eff) & (sum_exp > 0), 0.0)


def _scan_own_gene(aso_df, gene_to_data, cutoffs, n_jobs):
    """Scan each ASO against its own canonical gene."""
    return scan_gene_sites(
        aso_df=aso_df,
        gene_to_data=gene_to_data,
        target_genes=aso_df[CANONICAL_GENE_NAME].dropna().unique(),
        row_genes=aso_df[CANONICAL_GENE_NAME],
        cutoffs=cutoffs,
        n_jobs=n_jobs,
    )


def add_on_target_site_features(aso_df, gene_to_data, cutoffs, n_jobs=1):
    """Scores each ASO vs its own canonical gene -> BOTH on-target features per
    cutoff: ``on_target_total_hybridization_{c}`` and ``on_target_log_number_of_sites_{c}``.

    Returns (aso_df, [feature_names]).
    """
    sites = _scan_own_gene(aso_df, gene_to_data, cutoffs, n_jobs)
    columns = {}
    for cutoff in (int(c) for c in cutoffs):
        columns[f"on_target_total_hybridization_{cutoff}"] = total_hybridization(sites, cutoff)
        columns[f"on_target_log_number_of_sites_{cutoff}"] = log_number_of_sites(sites, cutoff)
    return add_columns(aso_df, columns), list(columns)


def on_target_total_hybridization(aso_df, gene_to_data, cutoffs, n_jobs=1):
    """Score each oligo against its own canonical gene: total hybridization, one feature per cutoff."""
    sites = _scan_own_gene(aso_df, gene_to_data, cutoffs, n_jobs)
    columns = {f"on_target_total_hybridization_{c}": total_hybridization(sites, c) for c in map(int, cutoffs)}
    return add_columns(aso_df, columns), list(columns)


def off_target_single_gene_hybridization(aso_df, gene_name, gene_to_data, cutoffs, n_jobs=1):
    """Score each oligo against a single fixed gene: total hybridization, one feature per cutoff."""
    sites = scan_gene_sites(
        aso_df=aso_df,
        gene_to_data=gene_to_data,
        target_genes=[gene_name],
        row_genes=[gene_name] * len(aso_df),
        cutoffs=cutoffs,
        n_jobs=n_jobs,
    )
    columns = {f"off_target_single_{gene_name}_c{c}": total_hybridization(sites, c) for c in map(int, cutoffs)}
    return add_columns(aso_df, columns), list(columns)
