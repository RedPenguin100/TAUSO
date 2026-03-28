import os
import uuid

import numpy as np
import pandas as pd

from ...data.consts import CANONICAL_GENE, SEQUENCE
from ...util import get_antisense
from ..hybridization.fast_hybridization import (
    TMP_PATH,
    Interaction,
    dump_target_file,
    get_trigger_mfe_scores_by_risearch,
)
from .off_target_functions import parse_risearch_output


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

    # 1. Setup Parallelization
    use_parallel = n_jobs != 1
    if use_parallel:
        try:
            from pandarallel import pandarallel

            verbose_score = 2 if verbose else 0
            pandarallel.initialize(
                nb_workers=n_jobs, progress_bar=verbose, verbose=verbose_score
            )
        except ImportError:
            if verbose:
                print("Warning: 'pandarallel' not found. Falling back to single core.")
            use_parallel = False

    RT = 0.616
    TMP_PATH.mkdir(exist_ok=True)

    # 2. Pre-compute target files for all required genes
    gene_to_target_info = {}

    try:
        for gene in target_genes:
            if gene in gene_to_data:
                sequence = gene_to_data[gene].full_mrna
                name_to_seq = {gene: sequence}

                target_hash = uuid.uuid4().hex
                target_filename = f"target-{gene}-{target_hash}.fa"
                target_path = dump_target_file(target_filename, name_to_seq)

                gene_to_target_info[gene] = {
                    "name_to_seq": name_to_seq,
                    "target_path": target_path,
                }
            else:
                if verbose:
                    print(f"Warning: Gene {gene} missing from gene_to_data.")

        # 3. Define the scoring function per row
        def calculate_score(row):
            gene = get_gene_fn(row)

            # Return 0.0 if the gene is NaN or wasn't successfully cached
            if pd.isna(gene) or gene not in gene_to_target_info:
                return 0.0

            aso_seq = row[SEQUENCE]
            trigger = get_antisense(aso_seq)
            run_id = f"{os.getpid()}-{uuid.uuid4().hex}"

            info = gene_to_target_info[gene]

            result_dict = get_trigger_mfe_scores_by_risearch(
                trigger,
                info["name_to_seq"],
                minimum_score=cutoff,
                parsing_type="2",
                interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
                transpose=True,
                target_file_cache=info["target_path"],
                unique_id=run_id,
            )

            result_df = parse_risearch_output(result_dict)

            if not result_df.empty and "energy" in result_df.columns:
                return np.exp(-RT * result_df["energy"]).sum()
            else:
                return 0.0

        # 4. Apply the function
        if use_parallel:
            scores = aso_df.parallel_apply(calculate_score, axis=1)
        else:
            scores = aso_df.apply(calculate_score, axis=1)

    finally:
        # 5. Cleanup all cached target FASTA files securely
        for info in gene_to_target_info.values():
            if os.path.exists(info["target_path"]):
                os.remove(info["target_path"])

    aso_df[feature_name] = scores
    return aso_df, feature_name


def on_target_total_hybridization(
        aso_df, gene_to_data, cutoff, n_jobs=1, verbose=False
):
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


def off_target_specific_seq_pandarallel(
        aso_df, gene_name, gene_to_data, cutoff, n_jobs=1, verbose=False
):
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
