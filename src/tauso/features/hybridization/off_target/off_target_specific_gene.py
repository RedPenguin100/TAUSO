"""Public per-gene hybridization features: score ASOs against a single target gene.

Thin entry points over the gene-chunk scoring orchestration — one for on-target (each
ASO vs its own canonical gene) and one for off-target against a fixed named gene. Each
returns the dataframe plus the names of the feature columns it added (one per cutoff).
"""

import math

import pandas as pd

from ....data.consts import CANONICAL_GENE_NAME
from ..fast_hybridization import sum_and_min_energy_by_trigger_multi_cutoff
from .gene_chunk_scoring import dispatch_gene_chunk_scores

_RT = 0.616


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


def _build_on_target_multiplicity_columns(aso_df, gene_scores, gene_to_row_triggers, cutoffs, RT=_RT):
    """From per-trigger ``(sum_exp, min_energy)``, emit ``on_target_total_hybridization_{c}`` (the Boltzmann
    sum, unchanged) and ``on_target_log_number_of_sites_{c}`` = ``log(sum_exp) - (-Emin/RT)`` (0 for one
    dominant site, log(N) for N equally strong sites)."""
    feature_names = []
    for cutoff in cutoffs:
        total = pd.Series(0.0, index=aso_df.index)
        logeff = pd.Series(0.0, index=aso_df.index)
        for gene, row_triggers in gene_to_row_triggers.items():
            gene_cutoff_scores = gene_scores[gene][cutoff]
            for idx, _ in row_triggers:
                value = gene_cutoff_scores.get(str(idx))
                if value is None:
                    continue
                sum_exp, min_energy = value
                total[idx] = sum_exp
                if sum_exp > 0:
                    le = math.log(max(sum_exp, 1e-300)) - (-min_energy / RT)
                    logeff[idx] = le if math.isfinite(le) else 0.0
        aso_df[f"on_target_total_hybridization_{cutoff}"] = total
        aso_df[f"on_target_log_number_of_sites_{cutoff}"] = logeff
        feature_names.extend([f"on_target_total_hybridization_{cutoff}", f"on_target_log_number_of_sites_{cutoff}"])
    return feature_names


def on_target_hybridization_with_multiplicity(aso_df, gene_to_data, cutoffs, n_jobs=1, RT=_RT):
    """On-target total hybridization AND site multiplicity, from ONE RIsearch pass per gene: emits
    ``on_target_total_hybridization_{c}`` (unchanged sum) and ``on_target_log_number_of_sites_{c}``."""
    return dispatch_gene_chunk_scores(
        aso_df=aso_df,
        gene_to_data=gene_to_data,
        target_genes=aso_df[CANONICAL_GENE_NAME].dropna().unique(),
        get_gene_fn=lambda row: row[CANONICAL_GENE_NAME],
        feature_name_fn=None,
        cutoffs=cutoffs,
        n_jobs=n_jobs,
        aggregation_factory=lambda cs: sum_and_min_energy_by_trigger_multi_cutoff(cs, RT=RT),
        build_columns=lambda df, gs, gtr, cs: _build_on_target_multiplicity_columns(df, gs, gtr, cs, RT),
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
