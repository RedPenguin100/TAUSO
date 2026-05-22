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
    """Pure logic function: Takes one row -> Returns one score."""
    trigger = get_antisense(row[SEQUENCE])
    target_gene = row[CANONICAL_GENE]

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

    result_df_agg = result_df_agg[result_df_agg["target"] != target_gene]
    if result_df_agg.empty:
        logger.warning(
            "Target gene %s excluded all off-target hits; score set to 0 for ASO %s",
            target_gene,
            row[SEQUENCE],
        )
        return 0

    simple_energies = dict(zip(result_df_agg["target"], result_df_agg["energy"]))
    return calculate_score_helper(simple_energies, general_exp_map, method)


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


def _parse_and_filter_hits(raw_output: str, group_df) -> dict:
    """Parse RIsearch TSV, min-aggregate per (trigger, target), remove self-hits.

    Returns {trigger_id: DataFrame of off-target hits}, empty dict if no valid hits.
    """
    all_hits = parse_risearch_output(raw_output)
    if all_hits.empty or "energy" not in all_hits.columns:
        return {}

    agg = all_hits.groupby(["trigger", "target"], sort=False)["energy"].min().reset_index()
    del all_hits

    idx_to_gene = {str(i): g for i, g in zip(group_df.index, group_df[CANONICAL_GENE])}
    agg["_own"] = agg["trigger"].map(idx_to_gene)
    filtered = agg[agg["target"] != agg["_own"]]
    del agg

    return {k: v for k, v in filtered.groupby("trigger", sort=False)}


def _score_triggers(hits_by_trigger: dict, indices, exp_map, method) -> pd.Series:
    """Apply calculate_score_helper for each trigger, return a Series indexed by indices."""
    scores = pd.Series(0.0, index=indices, dtype=float)
    for idx in indices:
        hits = hits_by_trigger.get(str(idx))
        if hits is not None:
            scores[idx] = calculate_score_helper(dict(zip(hits["target"], hits["energy"])), exp_map, method)
    return scores


def compute_group_batch(group_df, seq_map, exp_map, cutoff, method, prebuilt_target_path=None):
    """Batch version of compute_single_row: one RIsearch call for the entire group.

    prebuilt_target_path: if provided, the caller owns the file lifecycle.
    """
    if group_df.empty:
        return pd.Series(dtype=float)

    TMP_PATH.mkdir(exist_ok=True)

    _owns_target = prebuilt_target_path is None
    target_path = (
        dump_target_file(f"target-batch-{uuid.uuid4().hex}.fa", seq_map) if _owns_target else prebuilt_target_path
    )

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

    if not result.strip():
        return pd.Series(0.0, index=group_df.index, dtype=float)

    hits_by_trigger = _parse_and_filter_hits(result, group_df)
    del result

    return _score_triggers(hits_by_trigger, indices, exp_map, method)
