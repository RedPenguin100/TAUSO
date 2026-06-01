"""Off-target feature scores.

Calculates the various off-target feature scores from RIsearch hits, on top of the
abstracted RIsearch streaming utilities in fast_hybridization. Given a group of ASOs and
a prebuilt target FASTA, it runs RIsearch once, reduces the hits to the strongest energy
per (trigger, target) for each energy cutoff, and turns those energies into an
expression-weighted score.
"""

import logging
import os
import uuid

import pandas as pd

from ....data.consts import ASO_SEQUENCE, CANONICAL_GENE_NAME
from ....util import get_antisense
from ..fast_hybridization import (
    Interaction,
    min_energy_by_pair_multi_cutoff,
    parse_risearch_hits_pyarrow,
)

logger = logging.getLogger(__name__)


class AggregationMethod:
    ARTM = "ARTM"


def aggregate_energies(energy_dict, expression_dict, method):
    """Combine per-gene off-target energies into one expression-weighted score.

    energy_dict: {gene: energy} for the off-target genes a trigger binds.
    expression_dict: {gene: (TPM, normalised)}.
    method: an AggregationMethod. Only ARTM is supported: sum of energy * TPM-fraction
    over off-targets the trigger actually binds (energy < 0) and that are expressed.
    """
    if method != AggregationMethod.ARTM:
        raise ValueError(f"Unsupported aggregation method: {method}")
    if not energy_dict:
        return 0.0

    valid_targets = []
    for gene, energy in energy_dict.items():
        expr_tpm = expression_dict[gene][0] / 1e6
        if energy < 0 and expr_tpm > 0:
            valid_targets.append((gene, energy, expr_tpm))

    # Sort by gene name so summation order is identical regardless of dict/thread ordering.
    valid_targets.sort(key=lambda x: x[0])

    score = 0.0
    for _, energy, expr_tpm in valid_targets:
        score += energy * expr_tpm
    return score


def calculate_risearch_energy_per_cutoff(query_pairs, target_path, cutoffs, minimum_score):
    """Run RIsearch once at `minimum_score` and bucket the hits by cutoff.

    query_pairs: [(trigger_id, trigger_seq)]; target_path: a prebuilt target FASTA.
    cutoffs: the energy-score cutoffs to keep buckets for (a hit counts toward cutoff c
    when its score > c). minimum_score: the RIsearch -s threshold for the single run
    (the caller passes the loosest cutoff so every requested cutoff is derivable).

    Returns {cutoff: {(trigger, target): strongest_energy}} (no self-filter, no scoring).
    """
    if not query_pairs:
        return {int(c): {} for c in cutoffs}
    return parse_risearch_hits_pyarrow(
        trigger_id_seq_pairs=query_pairs,
        target_file_path=target_path,
        aggregation=min_energy_by_pair_multi_cutoff(cutoffs),
        minimum_score=minimum_score,
        parsing_type="2",
        interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
        transpose=True,
        batch_id=f"{os.getpid()}-{uuid.uuid4().hex}",
    )


def energy_score_per_cutoff(per_cutoff, indices, idx_to_gene, exp_map, method, cutoffs):
    """Turn per-cutoff off-target energies into a score Series per cutoff.

    per_cutoff: {cutoff: {(trigger, target): energy}} from calculate_risearch_energy_per_cutoff.
    indices: the ASO row indices to produce scores for (the Series index).
    idx_to_gene: {trigger_id: own canonical gene} — hits to a trigger's own gene are dropped.
    exp_map: {gene: (TPM, normalised)}; method: an AggregationMethod.

    Returns {cutoff: Series of scores indexed by `indices`}.
    """
    out = {}
    for cutoff in cutoffs:
        per_trigger: dict = {}
        for (trigger, target), energy in per_cutoff.get(int(cutoff), {}).items():
            if target == idx_to_gene.get(trigger):
                continue
            per_trigger.setdefault(trigger, {})[target] = energy

        scores = pd.Series(0.0, index=indices, dtype=float)
        for idx in indices:
            energies = per_trigger.get(str(idx))
            if energies:
                scores[idx] = aggregate_energies(energies, exp_map, method)
        out[cutoff] = scores
    return out


def compute_group_batch_multi_cutoff(group_df, exp_map, cutoffs, method, prebuilt_target_path):
    """Score one ASO group against a prebuilt target for every cutoff, with one RIsearch run.

    Returns {cutoff: Series of scores indexed like group_df}.
    """
    indices = group_df.index.tolist()
    if group_df.empty:
        return {c: pd.Series(dtype=float) for c in cutoffs}

    query_pairs = [(str(idx), get_antisense(seq)) for idx, seq in zip(indices, group_df[ASO_SEQUENCE])]
    idx_to_gene = {str(idx): gene for idx, gene in zip(indices, group_df[CANONICAL_GENE_NAME])}

    per_cutoff = calculate_risearch_energy_per_cutoff(query_pairs, prebuilt_target_path, cutoffs, min(cutoffs))
    return energy_score_per_cutoff(per_cutoff, indices, idx_to_gene, exp_map, method, cutoffs)
