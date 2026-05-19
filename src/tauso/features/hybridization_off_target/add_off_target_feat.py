import logging
import os
import uuid

import numpy as np
import pandas as pd

from ...data.consts import CANONICAL_GENE, SEQUENCE

logger = logging.getLogger(__name__)
from ...util import get_antisense
from ..hybridization.fast_hybridization import (
    TMP_PATH,
    Interaction,
    dump_target_file,
    get_trigger_mfe_scores_by_risearch,
    get_triggers_mfe_scores_batch,
)
from .off_target_functions import aggregate_off_targets, parse_risearch_output


class AggregationMethod:
    ARTM_log = "ARTM_log"
    ARTM = "ARTM"
    ARTM_weighted = "ARTM_weighted"
    GEO = "GEO"
    RANKED = "RANKED"
    MECH = "MECH"


def compute_single_row(row, general_seq_map, general_exp_map, cutoff, method):
    """
    Pure logic function: Takes one row -> Returns one score.
    """
    aso_seq = row[SEQUENCE]

    # 1. Get the target gene from the row
    target_gene = row[CANONICAL_GENE]

    trigger = get_antisense(aso_seq)

    # --- RIsearch ---
    result_dict = get_trigger_mfe_scores_by_risearch(
        trigger,
        general_seq_map,
        minimum_score=cutoff,
        interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
        parsing_type="2",
        transpose=True,
    )

    result_df = parse_risearch_output(result_dict)
    result_df_agg = aggregate_off_targets(result_df)

    if result_df_agg.empty:
        return 0

    # 2. Exclude the intended target from the results
    # We filter out any hit where the 'target' column matches the ASO's canonical gene
    result_df_agg = result_df_agg[result_df_agg["target"] != target_gene]

    if result_df_agg.empty:
        logger.warning(
            "Target gene %s excluded all off-target hits; score set to 0 for ASO %s", target_gene, row[SEQUENCE]
        )
        return 0

    # Extract energies and compute score
    # simple_energies: {Gene: Energy}
    simple_energies = dict(zip(result_df_agg["target"], result_df_agg["energy"]))

    score = calculate_score_helper(simple_energies, general_exp_map, method)
    return score


def calculate_score_helper(energy_dict, expression_dict, method):
    """
    Standardizes the scoring math to avoid code duplication.
    """
    score = 0.0
    RT = 0.616

    if not energy_dict:
        return 0.0

    # Pre-filter for valid interactions (Energy < 0) and Expression > 0
    # This speeds up the loops significantly
    valid_targets = []
    for gene, energy in energy_dict.items():
        expr_tpm = expression_dict[gene][0]
        # Use 1e6 scaling for TPM-based methods (1, 2, 3, 4, 5)
        # Method 0 uses norm, which shouldn't be scaled by 1e6
        if method != AggregationMethod.ARTM_log:
            expr_tpm = expr_tpm / 1e6

        if energy < 0 and expr_tpm > 0:
            valid_targets.append((gene, energy, expr_tpm))

    if not valid_targets:
        return 0.0

    if method == AggregationMethod.ARTM_log:  # Normalized Arithmetic (Direct exp, no 1e6 scaling)
        # Re-loop because we scaled exp above, but method 0 usually takes raw norm
        score = 0.0
        for gene, energy in energy_dict.items():
            expr_log = expression_dict.get(gene, 0)[1]  # Raw Norm
            if energy < 0 and expr_log > 0:
                score += energy * expr_log

    elif method == AggregationMethod.ARTM:  # TPM Arithmetic
        for _, energy, expr_tpm in valid_targets:
            score += energy * expr_tpm

    elif method == AggregationMethod.ARTM_weighted:  # Weighted Arithmetic (e^2)
        for _, energy, expr_tpm in valid_targets:
            score += energy * (expr_tpm**2)

    elif method == AggregationMethod.GEO:  # Geometric
        log_sum = 0.0
        for _, energy, expr_tpm in valid_targets:
            log_sum += expr_tpm * np.log(-energy)
        score = log_sum

    elif method == AggregationMethod.RANKED:  # Ranked
        # Sort by expression descending
        valid_targets.sort(key=lambda x: x[2], reverse=True)

        for rank, (_, energy, expr_tpm) in enumerate(valid_targets, start=1):
            batch = (rank - 1) // 10
            weight = 10 / (2**batch)
            score += energy * expr_tpm * weight

    elif method == AggregationMethod.MECH:  # Mechanical Statistics
        for _, energy, expr_tpm in valid_targets:
            score += expr_tpm * np.exp(-energy / RT)

    # Normalize Mechanical Statistics if required (usually done outside, but can be done here if context allows)
    # Note: Normalization requires knowing the max of the whole vector, so strictly speaking
    # method 5 normalization should happen after collecting all scores.
    # For now, this returns the raw Sum.

    return score


def compute_group_batch(group_df, seq_map, exp_map, cutoff, method, prebuilt_target_path=None):
    """Batch version of compute_single_row: one RIsearch call for the entire group.

    Pre-creates the target FASTA once from seq_map (or uses prebuilt_target_path), batches
    all queries, then applies the same per-row aggregation logic as compute_single_row.
    Results are identical to per-row scoring.

    prebuilt_target_path: if provided, the caller owns the file lifecycle (no create/delete here).
    """
    if group_df.empty:
        return pd.Series(dtype=float)

    TMP_PATH.mkdir(exist_ok=True)

    _owns_target = prebuilt_target_path is None
    target_path = (
        dump_target_file(f"target-batch-{uuid.uuid4().hex}.fa", seq_map) if _owns_target else prebuilt_target_path
    )

    # Avoid iterrows: build (id, antisense) pairs from column arrays
    indices = group_df.index.tolist()
    sequences = group_df[SEQUENCE].tolist()
    query_pairs = [(str(idx), get_antisense(seq)) for idx, seq in zip(indices, sequences)]

    try:
        result = get_triggers_mfe_scores_batch(
            trigger_id_seq_pairs=query_pairs,
            target_file_path=target_path,
            minimum_score=cutoff,
            parsing_type="2",
            interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
            transpose=True,
            batch_id=f"{os.getpid()}-{uuid.uuid4().hex}",
        )
    finally:
        if _owns_target and os.path.exists(target_path):
            os.remove(target_path)

    scores = pd.Series(0.0, index=group_df.index, dtype=float)

    if not result.strip():
        return scores

    all_hits = parse_risearch_output(result)
    del result

    if all_hits.empty or "energy" not in all_hits.columns:
        return scores

    # One aggregate pass: min energy per (trigger, target) — replaces N per-row groupby calls
    agg_all = all_hits.groupby(["trigger", "target"], sort=False)["energy"].min().reset_index()
    del all_hits

    # Vectorized self-target filter: map each trigger id → its canonical gene
    idx_to_gene = {str(i): g for i, g in zip(group_df.index, group_df[CANONICAL_GENE])}
    agg_all["_own"] = agg_all["trigger"].map(idx_to_gene)
    agg_filtered = agg_all[agg_all["target"] != agg_all["_own"]]
    del agg_all

    # O(1) lookup per trigger
    hits_by_trigger = {k: v for k, v in agg_filtered.groupby("trigger", sort=False)}
    del agg_filtered

    for idx in indices:
        trigger_id = str(idx)
        hits = hits_by_trigger.get(trigger_id)
        if hits is None:
            continue
        simple_energies = dict(zip(hits["target"], hits["energy"]))
        scores[idx] = calculate_score_helper(simple_energies, exp_map, method)

    return scores
