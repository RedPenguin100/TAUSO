"""Off-target feature scores.

Calculates the various off-target feature scores from RIsearch hits, on top of the
abstracted RIsearch streaming utilities in fast_hybridization. Given a group of ASOs and
a prebuilt target FASTA, it runs RIsearch once, reduces the hits to the strongest energy
per (trigger, target) for each energy cutoff, and turns those energies into an
expression-weighted score.
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
    PairAggregation,
    aggregate_by_pair_multi_cutoff,
    parse_risearch_hits_pyarrow,
)

logger = logging.getLogger(__name__)


class AggregationMethod:
    BOLTZMANN_SUM = "BOLTZ"


def aggregate_per_gene_tpm_weighted(energy_dict, expression_dict, method):
    """Combine per-gene off-target values into a log-abundance × log-free-energy score.

    energy_dict: {gene: Z_g} for the off-target genes a trigger binds; ``Z_g`` is the
    per-pair Boltzmann sum from PairAggregation.BOLTZMANN_SUM
    (Σ_sites exp(-energy/RT) for that (trigger, gene)).
    expression_dict: {gene: (TPM, normalised)}.
    method: AggregationMethod tag; consumed only by serialize_feature_name (kept in
    the signature so the existing call chain doesn't change).

    Computes Σ_g log1p(TPM_g) · log(Z_g). log(Z_g) is proportional to -ΔG_binding/RT
    for total off-target capture on transcript g; log1p(TPM_g) is a log-abundance
    weight (positive for any TPM > 0, sane magnitudes — housekeeping TPM~10⁴ gives
    ~9.2, low-expression TPM~10 gives ~2.4). Both quantities sit in a single-digit
    to ~50 range, so the per-gene contribution can't blow up the way TPM · Z would
    after the RT fix made Z values span ~10⁷ to ~10²¹.
    """
    if not energy_dict:
        return 0.0

    valid_targets = []
    for gene, energy in energy_dict.items():
        expr_tpm = expression_dict[gene][0]
        if expr_tpm > 0:
            valid_targets.append((gene, energy, expr_tpm))

    # Sort by gene name so summation order is identical regardless of dict/thread ordering.
    valid_targets.sort(key=lambda x: x[0])

    score = 0.0
    for _, energy, expr_tpm in valid_targets:
        score += math.log1p(expr_tpm) * math.log(energy)
    return score


def calculate_risearch_energy_per_cutoff(query_pairs, target_path, cutoffs, minimum_score):
    """Run RIsearch once at `minimum_score` and bucket the hits by cutoff.

    query_pairs: [(trigger_id, trigger_seq)]; target_path: a prebuilt target FASTA.
    cutoffs: the energy-score cutoffs to keep buckets for (a hit counts toward cutoff c
    when its score > c). minimum_score: the RIsearch -s threshold for the single run
    (the caller passes the loosest cutoff so every requested cutoff is derivable).

    Returns {cutoff: {(trigger, target): value}} where ``value`` is the per-pair
    reduction selected by PairAggregation.BOLTZMANN_SUM (no self-filter, no scoring).
    """
    if not query_pairs:
        return {int(c): {} for c in cutoffs}
    return parse_risearch_hits_pyarrow(
        trigger_id_seq_pairs=query_pairs,
        target_file_path=target_path,
        aggregation=aggregate_by_pair_multi_cutoff(cutoffs, strategy=PairAggregation.BOLTZMANN_SUM),
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
                scores[idx] = aggregate_per_gene_tpm_weighted(energies, exp_map, method)
        out[cutoff] = scores
    return out


def compute_group_batch_multi_cutoff_multi_topn(group_df, top_n_to_data, cutoffs, method, prebuilt_target_path):
    """Score one ASO group for several top_n levels from ONE RIsearch run.

    top_n_to_data: {top_n: (exp_map, gene_set)} where gene_set is the set of target
    genes that belong to head(top_n) of the cell-line / general expression table.
    Every gene_set must be a subset of the genes in ``prebuilt_target_path`` (which
    should be built at max(top_n)). Per-gene Boltzmann sums for smaller top_n are
    derived from the same RIsearch output by restricting to that top_n's gene_set —
    valid because RIsearch scores each (query, target) pair independently
    (see tests/complete/test_off_target_derivation_equivalence.py).

    Returns {(top_n, cutoff): Series of scores indexed like group_df}.
    """
    indices = group_df.index.tolist()
    if group_df.empty:
        return {(top_n, c): pd.Series(dtype=float) for top_n in top_n_to_data for c in cutoffs}

    query_pairs = [(str(idx), get_antisense(seq)) for idx, seq in zip(indices, group_df[ASO_SEQUENCE])]
    idx_to_gene = {str(idx): gene for idx, gene in zip(indices, group_df[CANONICAL_GENE_NAME])}

    per_cutoff = calculate_risearch_energy_per_cutoff(query_pairs, prebuilt_target_path, cutoffs, min(cutoffs))

    out: dict = {}
    for top_n, (exp_map, gene_set) in top_n_to_data.items():
        filtered = {
            c: {(trig, tgt): val for (trig, tgt), val in pairs.items() if tgt in gene_set}
            for c, pairs in per_cutoff.items()
        }
        for cutoff, series in energy_score_per_cutoff(filtered, indices, idx_to_gene, exp_map, method, cutoffs).items():
            out[(top_n, cutoff)] = series
    return out
