"""Off-target feature scores.

Runs RIsearch for a group of ASOs against a prebuilt target FASTA, reduces the hits to a
per-(trigger, target) Boltzmann occupancy score per cutoff, and turns those into an
expression-weighted off-target score.
"""

import logging
import math
import os
import uuid

import pandas as pd

from ....data.consts import ASO_SEQUENCE, CANONICAL_GENE_NAME
from ....util import get_antisense
from ..fast_hybridization import (
    Interaction,
    aggregate_by_pair_multi_cutoff,
    parse_risearch_hits_pyarrow,
)

logger = logging.getLogger(__name__)


class AggregationMethod:
    BOLTZMANN_SUM = "BOLTZ"


def aggregate_per_gene_tpm_weighted(occupancy_score_by_gene, expression_map, method):
    """Combine per-gene off-target occupancy scores into one expression-weighted score.

    occupancy_score_by_gene: {gene: Z_g}, the per-(trigger, gene) Boltzmann sum Σ_sites exp(-energy/RT).
    expression_map: {gene: (TPM, normalised)}. method is unused (kept for the shared signature).

    Returns Σ_g log1p(TPM_g) · log(Z_g): log(Z_g) tracks -ΔG_binding/RT for total capture on gene
    g, log1p(TPM_g) is a log-abundance weight; both stay small so no single gene dominates.
    """
    if not occupancy_score_by_gene:
        return 0.0

    valid_targets = []
    for gene, occupancy_score in occupancy_score_by_gene.items():
        expression_tpm = expression_map[gene][0]
        if expression_tpm > 0:
            valid_targets.append((gene, occupancy_score, expression_tpm))

    # Sort by gene name so the summation order is independent of dict/thread ordering.
    valid_targets.sort(key=lambda item: item[0])

    score = 0.0
    for _, occupancy_score, expression_tpm in valid_targets:
        score += math.log1p(expression_tpm) * math.log(occupancy_score)
    return score


def risearch_occupancy_score_per_cutoff(query_pairs, target_path, cutoffs, minimum_score):
    """Run RIsearch once and reduce the hits to a per-(trigger, target) Boltzmann occupancy score.

    query_pairs: [(trigger_id, trigger_seq)]; target_path: a prebuilt target FASTA. minimum_score
    is the RIsearch -s threshold (the loosest cutoff, so every cutoff in `cutoffs` is derivable —
    a hit counts toward a cutoff when its score exceeds it).

    Returns {cutoff: {(trigger, target): Z}}, Z = Σ_sites exp(-energy/RT) (no self-filter).
    """
    if not query_pairs:
        return {int(cutoff): {} for cutoff in cutoffs}
    return parse_risearch_hits_pyarrow(
        trigger_id_seq_pairs=query_pairs,
        target_file_path=target_path,
        aggregation=aggregate_by_pair_multi_cutoff(cutoffs),
        minimum_score=minimum_score,
        parsing_type="2",
        interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
        transpose=True,
        batch_id=f"{os.getpid()}-{uuid.uuid4().hex}",
    )


def offtarget_score_per_cutoff(per_cutoff, aso_indices, gene_by_trigger, expression_map, method, cutoffs):
    """Turn per-cutoff off-target occupancy scores into a score Series per cutoff.

    per_cutoff: {cutoff: {(trigger, target): Z}} from risearch_occupancy_score_per_cutoff.
    gene_by_trigger: {trigger_id: own canonical gene} — hits to a trigger's own gene are dropped.
    expression_map: {gene: (TPM, normalised)}; method is passed through unused.

    Returns {cutoff: Series of scores indexed by `aso_indices`}.
    """
    scores_by_cutoff = {}
    for cutoff in cutoffs:
        per_trigger: dict = {}
        for (trigger, target), occupancy_score in per_cutoff.get(int(cutoff), {}).items():
            if target == gene_by_trigger.get(trigger):
                continue
            per_trigger.setdefault(trigger, {})[target] = occupancy_score

        scores = pd.Series(0.0, index=aso_indices, dtype=float)
        for aso_index in aso_indices:
            occupancy_scores = per_trigger.get(str(aso_index))
            if occupancy_scores:
                scores[aso_index] = aggregate_per_gene_tpm_weighted(occupancy_scores, expression_map, method)
        scores_by_cutoff[cutoff] = scores
    return scores_by_cutoff


def compute_group_batch_multi_cutoff_multi_topn(group_df, top_n_to_data, cutoffs, method, prebuilt_target_path):
    """Score one ASO group for several top_n levels from a single RIsearch pass.

    top_n_to_data: {top_n: (expression_map, gene_set)}; each gene_set is head(top_n) of the
    expression table and must be a subset of the genes in ``prebuilt_target_path`` (built at
    max(top_n)). Smaller top_n are derived by restricting the shared occupancy scores to that
    gene_set — valid because RIsearch scores each (query, target) pair independently.

    Returns {(top_n, cutoff): Series of scores indexed like group_df}.
    """
    aso_indices = group_df.index.tolist()
    if group_df.empty:
        return {(top_n, cutoff): pd.Series(dtype=float) for top_n in top_n_to_data for cutoff in cutoffs}

    query_pairs = [
        (str(aso_index), get_antisense(sequence)) for aso_index, sequence in zip(aso_indices, group_df[ASO_SEQUENCE])
    ]
    gene_by_trigger = {str(aso_index): gene for aso_index, gene in zip(aso_indices, group_df[CANONICAL_GENE_NAME])}

    per_cutoff = risearch_occupancy_score_per_cutoff(query_pairs, prebuilt_target_path, cutoffs, min(cutoffs))

    scores_by_top_n_cutoff: dict = {}
    for top_n, (expression_map, gene_set) in top_n_to_data.items():
        occupancy_score_in_top_n = {
            cutoff: {
                (trigger, target): occupancy_score
                for (trigger, target), occupancy_score in pairs.items()
                if target in gene_set
            }
            for cutoff, pairs in per_cutoff.items()
        }
        for cutoff, series in offtarget_score_per_cutoff(
            occupancy_score_in_top_n, aso_indices, gene_by_trigger, expression_map, method, cutoffs
        ).items():
            scores_by_top_n_cutoff[(top_n, cutoff)] = series
    return scores_by_top_n_cutoff
