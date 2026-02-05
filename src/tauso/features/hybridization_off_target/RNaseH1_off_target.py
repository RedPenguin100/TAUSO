import os
import uuid

import pandas as pd
import numpy as np

from .off_target_functions import parse_risearch_output
from ..hybridization.fast_hybridization import get_trigger_mfe_scores_by_risearch, Interaction, dump_target_file, \
    TMP_PATH
from ...new_model.consts_dataframe import SEQUENCE, CANONICAL_GENE
from ...util import get_antisense


def off_target_specific_seq_pandarallel(aso_df, gene_name, gene_to_data, cutoff, n_jobs=1, verbose=False):
    # 2. Setup Parallelization
    use_parallel = n_jobs != 1
    if use_parallel:
        try:
            from pandarallel import pandarallel
            verbose_score = 2 if verbose else 0
            pandarallel.initialize(nb_workers=n_jobs, progress_bar=verbose, verbose=verbose_score)
        except ImportError:
            if verbose: print("Warning: 'pandarallel' not found. Falling back to single core.")
            use_parallel = False

    sequence = gene_to_data[gene_name].full_mrna
    name_to_seq = {gene_name: sequence}
    RT = 0.616
    feature_name = f"off_target_single_{gene_name}_c{cutoff}"

    TMP_PATH.mkdir(exist_ok=True)
    target_hash = uuid.uuid4().hex
    target_filename = f"target-{target_hash}.fa"
    target_path = dump_target_file(target_filename, name_to_seq)

    def calculate_score(row):
        aso_seq = row[SEQUENCE]
        trigger = get_antisense(aso_seq)

        result_dict = get_trigger_mfe_scores_by_risearch(
            trigger,
            name_to_seq,
            minimum_score=cutoff,
            parsing_type='2',
            interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
            transpose=True,
            target_file_cache=target_path
        )

        result_df = parse_risearch_output(result_dict)

        if not result_df.empty and 'energy' in result_df.columns:
            return np.exp(-RT * result_df['energy']).sum()
        else:
            return 0.0

    try:
        if use_parallel:
            scores = aso_df.parallel_apply(calculate_score, axis=1)
        else:
            scores = aso_df.apply(calculate_score, axis=1)
    finally:
        if os.path.exists(target_path):
            os.remove(target_path)

    aso_df[feature_name] = scores

    return aso_df, feature_name


def calculate_selectivity_parallel(aso_df, off_target_gene_name, gene_to_data, cutoff=500, n_jobs=1, verbose=False):
    """
    Calculates the Selectivity Ratio (On-Target / (On-Target + Off-Target)) in one pass.

    Args:
        aso_df: DataFrame containing 'SEQUENCE' and 'CANONICAL_GENE'.
        off_target_gene_name: The specific gene (e.g., pseudogene) to test as the off-target.
        gene_to_data: Dictionary mapping gene names to objects with .full_mrna.

    Returns:
        aso_df (updated), feature_name (str)
    """

    # 1. Setup Parallelization
    use_parallel = n_jobs != 1
    if use_parallel:
        try:
            from pandarallel import pandarallel
            pandarallel.initialize(nb_workers=n_jobs, progress_bar=verbose, verbose=1 if verbose else 0)
        except ImportError:
            if verbose: print("Warning: 'pandarallel' not found. Falling back to single core.")
            use_parallel = False

    RT = 0.616
    TMP_PATH.mkdir(exist_ok=True)

    # --- OPTIMIZATION: PRE-WRITE TARGET FILES ---

    # A. Prepare the OFF-TARGET file (Fixed for this entire run)
    off_target_seq = gene_to_data[off_target_gene_name].full_mrna
    off_target_hash = uuid.uuid4().hex
    off_target_filename = f"target_OFF_{off_target_gene_name}-{off_target_hash}.fa"
    off_target_path = dump_target_file(off_target_filename, {off_target_gene_name: off_target_seq})

    # B. Prepare ON-TARGET files (Variable based on row)
    unique_on_genes = aso_df[CANONICAL_GENE].unique()
    on_target_paths = {}
    created_files = [off_target_path]

    if verbose:
        print(f"Pre-generating target files for {len(unique_on_genes)} unique on-target genes...")

    for gene in unique_on_genes:
        try:
            seq = gene_to_data[gene].full_mrna
            f_hash = uuid.uuid4().hex
            f_name = f"target_ON_{gene}-{f_hash}.fa"
            f_path = dump_target_file(f_name, {gene: seq})

            on_target_paths[gene] = f_path
            created_files.append(f_path)
        except KeyError:
            if verbose: print(f"Warning: Gene {gene} not found in gene_to_data. Rows will score 0.")
            on_target_paths[gene] = None

    # Define the Single Feature Name
    feature_name = f"selectivity_{off_target_gene_name}_c{cutoff}"

    # --- WORKER FUNCTION ---
    def process_row(row):
        aso_seq = row[SEQUENCE]
        on_gene = row[CANONICAL_GENE]

        trigger = get_antisense(aso_seq)

        # 1. Calculate ON-TARGET Score
        score_on = 0.0
        on_path = on_target_paths.get(on_gene)

        if on_path:
            res_on = get_trigger_mfe_scores_by_risearch(
                trigger, {on_gene: "dummy"},
                minimum_score=cutoff,
                parsing_type='2', interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
                transpose=True, target_file_cache=on_path
            )
            df_on = parse_risearch_output(res_on)
            if not df_on.empty and 'energy' in df_on.columns:
                score_on = np.exp(-RT * df_on['energy']).sum()

        # 2. Calculate OFF-TARGET Score
        score_off = 0.0
        res_off = get_trigger_mfe_scores_by_risearch(
            trigger, {off_target_gene_name: "dummy"},
            minimum_score=cutoff,
            parsing_type='2', interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
            transpose=True, target_file_cache=off_target_path
        )
        df_off = parse_risearch_output(res_off)
        if not df_off.empty and 'energy' in df_off.columns:
            score_off = np.exp(-RT * df_off['energy']).sum()

        # 3. Calculate Ratio
        # On / (On + Off) -> Closer to 1.0 is better
        return score_on / (score_on + score_off + 1e-9)

    # --- EXECUTION ---
    try:
        if use_parallel:
            scores = aso_df.parallel_apply(process_row, axis=1)
        else:
            scores = aso_df.apply(process_row, axis=1)

        aso_df[feature_name] = scores

    finally:
        # --- CLEANUP ---
        for f_path in created_files:
            if f_path and os.path.exists(f_path):
                os.remove(f_path)

    return aso_df, feature_name