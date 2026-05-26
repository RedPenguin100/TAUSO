"""Off-target reductions + scoring.

Personality: turn a group of ASOs into off-target feature scores. It runs one RIsearch
pass per group (via the generic streaming primitives in fast_hybridization), reduces the
hits to min-energy-per-(trigger, target) for each cutoff, then applies the off-target
scoring math (self-target filtering + the AggregationMethod weighting). No RIsearch
process / parsing mechanism lives here — that's fast_hybridization's job.
"""

import logging
import os
import uuid

import numpy as np
import pandas as pd

from ....data.consts import CANONICAL_GENE, SEQUENCE
from ....util import get_antisense
from ..fast_hybridization import (
    TMP_PATH,
    Interaction,
    min_energy_by_pair_multi_cutoff,
    parse_risearch_hits_pyarrow,
)

logger = logging.getLogger(__name__)


class AggregationMethod:
    ARTM_log = "ARTM_log"
    ARTM = "ARTM"
    ARTM_weighted = "ARTM_weighted"
    GEO = "GEO"
    RANKED = "RANKED"
    MECH = "MECH"


def calculate_score_helper(energy_dict, expression_dict, method):
    """Standardizes the scoring math to avoid code duplication."""
    if not energy_dict:
        return 0.0

    RT = 0.616
    score = 0.0

    valid_targets = []
    for gene, energy in energy_dict.items():
        expr_tpm = expression_dict[gene][0]
        if method != AggregationMethod.ARTM_log:
            expr_tpm = expr_tpm / 1e6
        if energy < 0 and expr_tpm > 0:
            valid_targets.append((gene, energy, expr_tpm))

    if not valid_targets:
        return 0.0

    # Sort by gene name so summation order is identical regardless of dict/thread ordering.
    valid_targets.sort(key=lambda x: x[0])

    if method == AggregationMethod.ARTM_log:
        for gene, energy in sorted(energy_dict.items()):
            expr_log = expression_dict.get(gene, 0)[1]
            if energy < 0 and expr_log > 0:
                score += energy * expr_log

    elif method == AggregationMethod.ARTM:
        for _, energy, expr_tpm in valid_targets:
            score += energy * expr_tpm

    elif method == AggregationMethod.ARTM_weighted:
        for _, energy, expr_tpm in valid_targets:
            score += energy * (expr_tpm**2)

    elif method == AggregationMethod.GEO:
        for _, energy, expr_tpm in valid_targets:
            score += expr_tpm * np.log(-energy)

    elif method == AggregationMethod.RANKED:
        valid_targets.sort(key=lambda x: x[2], reverse=True)
        for rank, (_, energy, expr_tpm) in enumerate(valid_targets, start=1):
            batch = (rank - 1) // 10
            weight = 10 / (2**batch)
            score += energy * expr_tpm * weight

    elif method == AggregationMethod.MECH:
        for _, energy, expr_tpm in valid_targets:
            score += expr_tpm * np.exp(-energy / RT)

    return score


def risearch_min_energy_multi_cutoff(query_pairs, target_path, cutoffs):
    """ONE RIsearch pass at the loosest cutoff against a prebuilt target.

    Returns {cutoff: {(trigger, target): min_energy}} (no self-filter, no scoring)."""
    if not query_pairs:
        return {int(c): {} for c in cutoffs}
    return parse_risearch_hits_pyarrow(
        trigger_id_seq_pairs=query_pairs,
        target_file_path=target_path,
        aggregation=min_energy_by_pair_multi_cutoff(cutoffs),
        minimum_score=min(cutoffs),
        parsing_type="2",
        interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
        transpose=True,
        batch_id=f"{os.getpid()}-{uuid.uuid4().hex}",
    )


def score_min_energy_per_cutoff(per_cutoff, indices, idx_to_gene, exp_map, method, cutoffs):
    """Self-filter + calculate_score_helper for each cutoff. Returns {cutoff: Series}."""
    out = {}
    for cutoff in cutoffs:
        min_energy = per_cutoff.get(int(cutoff), {})
        # Group by trigger, filter out self-targets, build {trigger: {target: energy}}
        trigger_to_energies: dict = {}
        for (trig, tgt), e in min_energy.items():
            if tgt == idx_to_gene.get(trig):
                continue  # self-hit
            trigger_to_energies.setdefault(trig, {})[tgt] = e

        scores = pd.Series(0.0, index=indices, dtype=float)
        for idx in indices:
            energies = trigger_to_energies.get(str(idx))
            if energies:
                scores[idx] = calculate_score_helper(energies, exp_map, method)
        out[cutoff] = scores
    return out


def compute_group_batch_multi_cutoff(group_df, exp_map, cutoffs, method, prebuilt_target_path):
    """Score a group against a prebuilt target for several cutoffs with ONE RIsearch pass.

    Runs RIsearch once at the loosest cutoff and derives every cutoff's min-energy-per-pair
    in the same streaming pass (see min_energy_by_pair_multi_cutoff and
    tests/complete/test_off_target_derivation_equivalence.py).

    Returns {cutoff: Series of scores indexed like group_df}.
    """
    indices = group_df.index.tolist()
    if group_df.empty:
        return {c: pd.Series(dtype=float) for c in cutoffs}

    TMP_PATH.mkdir(exist_ok=True)
    sequences = group_df[SEQUENCE].tolist()
    query_pairs = [(str(idx), get_antisense(seq)) for idx, seq in zip(indices, sequences)]
    idx_to_gene = {str(i): g for i, g in zip(group_df.index, group_df[CANONICAL_GENE])}

    per_cutoff = risearch_min_energy_multi_cutoff(query_pairs, prebuilt_target_path, cutoffs)
    return score_min_energy_per_cutoff(per_cutoff, indices, idx_to_gene, exp_map, method, cutoffs)
