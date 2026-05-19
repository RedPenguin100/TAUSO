import logging
import os
import uuid

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

from ...data.consts import CANONICAL_GENE, SEQUENCE
from ...util import get_antisense
from ..hybridization.fast_hybridization import (
    TMP_PATH,
    Interaction,
    dump_target_file,
    get_triggers_mfe_scores_batch,
)
from .off_target_functions import parse_risearch_output

_RT = 0.616


def _sum_exp_energy_by_trigger(result_df: pd.DataFrame) -> dict:
    """Return {trigger_id: sum(exp(-RT * energy))} for every trigger in result_df."""
    exp_vals = np.exp(-_RT * result_df["energy"].to_numpy())
    return result_df.assign(_exp=exp_vals).groupby("trigger", sort=False)["_exp"].sum().to_dict()


def _validate_genes_found(target_genes, gene_to_data):
    not_found_genes = []
    for gene in target_genes:
        if gene not in gene_to_data:
            not_found_genes.append(gene)
    if not_found_genes:
        raise ValueError(f"The following genes are not found in gene_to_data: {not_found_genes}")


def _apply_risearch_scoring(
    aso_df,
    gene_to_data,
    target_genes,
    get_gene_fn,
    feature_name,
    cutoff,
    n_jobs,
    verbose,
):
    """Core logic for RIsearch hybridization scoring to avoid code duplication."""
    _validate_genes_found(target_genes, gene_to_data)

    TMP_PATH.mkdir(exist_ok=True)

    # Pre-compute target files for all required genes (one per gene, reused for all rows)
    gene_to_target_info = {}

    try:
        for gene in target_genes:
            sequence = gene_to_data[gene].full_mrna
            target_path = dump_target_file(f"target-{gene}-{uuid.uuid4().hex}.fa", {gene: sequence})
            gene_to_target_info[gene] = {"target_path": target_path}

        # Group rows by their target gene so we fire one RIsearch call per gene
        # instead of one per row. Use column arrays to avoid per-row Series creation.
        indices = aso_df.index.tolist()
        seqs = aso_df[SEQUENCE].tolist()
        genes_col = aso_df.apply(get_gene_fn, axis=1).tolist()

        from collections import defaultdict

        gene_to_row_triggers = defaultdict(list)
        for idx, seq, gene in zip(indices, seqs, genes_col):
            if pd.isna(gene) or gene not in gene_to_target_info:
                continue
            gene_to_row_triggers[gene].append((idx, get_antisense(seq)))

        scores = pd.Series(0.0, index=aso_df.index)

        for gene, row_triggers in gene_to_row_triggers.items():
            result = get_triggers_mfe_scores_batch(
                trigger_id_seq_pairs=[(str(idx), trig) for idx, trig in row_triggers],
                target_file_path=gene_to_target_info[gene]["target_path"],
                minimum_score=cutoff,
                parsing_type="2",
                interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
                transpose=True,
                batch_id=f"{os.getpid()}-{uuid.uuid4().hex}",
            )

            if not result.strip():
                continue

            result_df = parse_risearch_output(result)
            del result

            if result_df.empty or "energy" not in result_df.columns:
                continue

            score_dict = _sum_exp_energy_by_trigger(result_df)
            del result_df
            for idx, _ in row_triggers:
                val = score_dict.get(str(idx))
                if val is not None:
                    scores[idx] = val

    finally:
        for info in gene_to_target_info.values():
            if os.path.exists(info["target_path"]):
                os.remove(info["target_path"])

    aso_df[feature_name] = scores
    return aso_df, feature_name


def on_target_total_hybridization(aso_df, gene_to_data, cutoff, n_jobs=1, verbose=False):
    """Scores against the dynamic canonical gene found in each row."""
    unique_genes = aso_df[CANONICAL_GENE].dropna().unique()
    feature_name = f"on_target_total_hybridization_{cutoff}"

    return _apply_risearch_scoring(
        aso_df=aso_df,
        gene_to_data=gene_to_data,
        target_genes=unique_genes,
        get_gene_fn=lambda row: row[CANONICAL_GENE],  # Fetch gene dynamically
        feature_name=feature_name,
        cutoff=cutoff,
        n_jobs=n_jobs,
        verbose=verbose,
    )


def off_target_specific_seq_pandarallel(aso_df, gene_name, gene_to_data, cutoff, n_jobs=1, verbose=False):
    """Scores against a single statically provided gene."""
    feature_name = f"off_target_single_{gene_name}_c{cutoff}"

    return _apply_risearch_scoring(
        aso_df=aso_df,
        gene_to_data=gene_to_data,
        target_genes=[gene_name],
        get_gene_fn=lambda row: gene_name,  # Fetch gene statically
        feature_name=feature_name,
        cutoff=cutoff,
        n_jobs=n_jobs,
        verbose=verbose,
    )
