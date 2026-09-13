"""Off-target feature scores.

Runs RIsearch for a group of ASOs against a prebuilt target FASTA, reduces the hits to a
per-(query, target) Boltzmann occupancy score per cutoff, and turns those into an
expression-weighted off-target score.
"""

import logging

import numpy as np
import pandas as pd
import pyrisearch_tauso

from ....data.consts import ASO_SEQUENCE, CANONICAL_GENE_NAME
from . import ASO_TARGET_MATRIX, EXTENSION_PENALTY, RT_KCAL_MOL

logger = logging.getLogger(__name__)


class AggregationMethod:
    BOLTZMANN_SUM = "BOLTZ"


def occupancy_hits(query_pairs, target_path, cutoffs):
    """Run RIsearch once and tabulate what each query hit, per cutoff.

    query_pairs: [(query_id, query_seq)]; target_path: a prebuilt target FASTA.
    The loosest cutoff is the search threshold; each cutoff keeps score > cutoff.

    One row per (cutoff, query, target) with the Boltzmann occupancy summed over
    that pair's sites. No self-filter.
    """
    per_cutoff = pyrisearch_tauso.search_reduced(
        queries=query_pairs,
        targets=target_path,
        reduction=pyrisearch_tauso.energy_stats(cutoffs, rt=RT_KCAL_MOL),
        min_score=min(cutoffs),
        matrix=ASO_TARGET_MATRIX,
        extension_penalty=EXTENSION_PENALTY,
        transpose=True,
    )
    return pd.DataFrame(
        [
            (cutoff, query, target, stats.sum_exp)
            for cutoff, pairs in per_cutoff.items()
            for (query, target), stats in pairs.items()
        ],
        columns=["cutoff", "query", "target", "occupancy"],
    )


def tpm_weighted_score(hits, expression_map):
    """Score each query: summed over its genes, log1p(TPM) * log(occupancy).

    log(occupancy) tracks -dG_binding/RT for total capture on a gene and
    log1p(TPM) is a log-abundance weight; both stay small, so no one gene
    dominates. A gene with no expression carries no weight and is left out.

    Sorted by gene so the total does not depend on the order the hits arrived
    in. Returns a Series indexed by query id; a query with nothing left is
    absent from it.
    """
    scored = hits.assign(tpm=hits["target"].map(expression_map))
    scored = scored[scored["tpm"] > 0].sort_values("target")
    if scored.empty:
        return pd.Series(dtype=float)
    terms = np.log1p(scored["tpm"]) * np.log(scored["occupancy"])
    return terms.groupby(scored["query"]).sum()


def compute_group_batch_multi_cutoff_multi_topn(group_df, top_n_to_data, cutoffs, prebuilt_target_path):
    """Score one ASO group for several top_n levels from a single RIsearch pass.

    top_n_to_data: {top_n: (expression_map, gene_set)}; each gene_set is head(top_n)
    of the expression table and must be a subset of the genes in
    ``prebuilt_target_path`` (built at max(top_n)). A smaller top_n keeps the hits
    whose target is in its gene_set, which is the same answer as searching that
    universe, because RIsearch scores each (query, target) pair independently.

    Returns {(top_n, cutoff): Series of scores indexed like group_df}.
    """
    cutoffs = [int(c) for c in cutoffs]
    if group_df.empty:
        return {(top_n, cutoff): pd.Series(dtype=float) for top_n in top_n_to_data for cutoff in cutoffs}

    query_of_row = {row: str(row) for row in group_df.index}
    hits = occupancy_hits(
        [(query_of_row[row], seq) for row, seq in zip(group_df.index, group_df[ASO_SEQUENCE])],
        prebuilt_target_path,
        cutoffs,
    )

    # An ASO hitting the gene it was designed against is on target, not off.
    own_gene = {query_of_row[row]: gene for row, gene in zip(group_df.index, group_df[CANONICAL_GENE_NAME])}
    if not hits.empty:
        hits = hits[hits["target"] != hits["query"].map(own_gene)]

    scores = {}
    for top_n, (expression_map, gene_set) in top_n_to_data.items():
        in_universe = hits[hits["target"].isin(gene_set)] if not hits.empty else hits
        for cutoff in cutoffs:
            by_query = tpm_weighted_score(in_universe[in_universe["cutoff"] == cutoff], expression_map)
            scores[top_n, cutoff] = pd.Series(
                [by_query.get(query_of_row[row], 0.0) for row in group_df.index],
                index=group_df.index,
                dtype=float,
            )
    return scores
