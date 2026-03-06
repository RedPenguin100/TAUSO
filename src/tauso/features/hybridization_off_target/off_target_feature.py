import numpy as np
import pandas as pd


from .add_off_target_feat import compute_single_row
from ...data.consts import CELL_LINE_DEPMAP


def serialize_feature_name(method, top_n, cutoff, is_specific):
    if is_specific:
        return f"off_target_score_specific_{method}_n{top_n}_c{cutoff}"
    return f"off_target_score_general_{method}_n{top_n}_c{cutoff}"


def populate_off_target_specific(
        ASO_df,
        gene_to_data,
        cell_line2data,
        top_n_list,
        cutoff_list,
        method,
        n_cores=1
):
    """
    Enriches ASO_df with off-target scores based on the specific cell line transcriptome.
    Uses CELL_LINE_DEPMAP column to map rows to cell_line2data keys.
    """
    ASO_df = ASO_df.copy()
    feature_names = []

    # Initialize pandarallel if needed
    if n_cores > 1:
        from pandarallel import pandarallel
        try:
            pandarallel.initialize(nb_workers=n_cores, verbose=0)
        except RuntimeError:
            pass  # Already initialized

    for top_n in top_n_list:
        for cutoff in cutoff_list:
            print(f"Running Specific Config: TopN={top_n}, Cutoff={cutoff}")
            spec_col_name = serialize_feature_name(method, top_n, cutoff, is_specific=True)

            # Container for the results
            # We will compute chunks and concat them back together
            calculated_scores = []

            # Group by Cell Line to apply specific transcriptomes
            grouped = ASO_df.groupby(CELL_LINE_DEPMAP)

            for cell_line, group_df in grouped:
                # 1. Validation: Check if we have data for this cell line
                if cell_line not in cell_line2data:
                    # If missing context, score is NaN
                    empty_res = pd.Series([np.nan] * len(group_df), index=group_df.index)
                    calculated_scores.append(empty_res)
                    continue

                # 2. Get Specific Transcriptome (Top N)
                specific_df = cell_line2data[cell_line].head(top_n)

                # 3. Build Maps for this specific cell line
                spec_seq_map = {}
                spec_exp_map = {}

                norm_col = next((c for c in specific_df.columns if 'expression_norm' in c), 'expression_norm')

                for _, row in specific_df.iterrows():
                    gene = row['Gene'].split()[0]
                    if gene in gene_to_data:
                        spec_seq_map[gene] = gene_to_data[gene].full_mrna
                        spec_exp_map[gene] = (
                            row.get('expression_TPM', 0),
                            row.get(norm_col, row.get('expression_norm', 0))
                        )

                # 4. Compute Scores for this Group
                # We assume the worker function is available in scope
                # Optimization: Only parallelize if group is large enough to justify overhead
                if n_cores > 1 and len(group_df) > n_cores:
                    group_scores = group_df.parallel_apply(
                        compute_single_row,  # Reusing the generic worker
                        axis=1,
                        args=(spec_seq_map, spec_exp_map, cutoff, method)
                    )
                else:
                    group_scores = group_df.apply(
                        compute_single_row,
                        axis=1,
                        args=(spec_seq_map, spec_exp_map, cutoff, method)
                    )

                calculated_scores.append(group_scores)

            # 5. Reassemble and Assign
            if calculated_scores:
                full_series = pd.concat(calculated_scores)
                # Map back to original dataframe indices to ensure alignment
                ASO_df[spec_col_name] = full_series.reindex(ASO_df.index)
            else:
                ASO_df[spec_col_name] = np.nan

            feature_names.append(spec_col_name)

    return ASO_df, feature_names


def populate_off_target_general(
        ASO_df,
        gene_to_data,
        cell_line2data,
        top_n_list,
        cutoff_list,
        method,
        n_cores=1  # Added parameter
):
    """
    Enriches ASO_df with off-target scores.
    Supports parallelization via pandarallel.
    """
    ASO_df = ASO_df.copy()
    feature_names = []

    if 'general' not in cell_line2data:
        raise ValueError("Key 'general' not found in cell_line2data dictionary.")
    general_df_all = cell_line2data['general']

    # Initialize pandarallel if needed
    if n_cores > 1:
        from pandarallel import pandarallel
        # verbose=1 shows a progress bar, use 0 to silence
        pandarallel.initialize(nb_workers=n_cores, verbose=1, progress_bar=True)

    for top_n in top_n_list:
        general_df = general_df_all.head(top_n)

        # 1. Prepare Maps ONCE
        general_seq_map = {}
        general_exp_map = {}

        norm_col = next((c for c in general_df.columns if 'expression_norm' in c), 'expression_norm')
        for _, row in general_df.iterrows():
            gene = row['Gene'].split()[0]
            if gene not in gene_to_data:
                raise KeyError(f"CRITICAL: Gene '{gene}' not found in gene_to_data mapping.")

            general_seq_map[gene] = gene_to_data[gene].full_mrna
            general_exp_map[gene] = (
                row.get('expression_TPM', 0),
                row.get(norm_col, row.get('expression_norm', 0))
            )

        for cutoff in cutoff_list:
            print(f"Running config: TopN={top_n}, Cutoff={cutoff}")
            gen_col_name = serialize_feature_name(method, top_n, cutoff, is_specific=False)

            # 2. Parallel vs Single-Core Execution
            if n_cores > 1:
                scores = ASO_df.parallel_apply(
                    compute_single_row,
                    axis=1,
                    args=(general_seq_map, general_exp_map, cutoff, method)
                )
            else:
                scores = ASO_df.apply(
                    compute_single_row,
                    axis=1,
                    args=(general_seq_map, general_exp_map, cutoff, method)
                )

            # Assign new column
            ASO_df[gen_col_name] = scores
            feature_names.append(gen_col_name)

    return ASO_df, feature_names


def populate_off_target_specific_per_rank(
        ASO_df,
        gene_to_data,
        cell_line2data,
        max_rank,
        cutoff_list,
        method,
        n_cores=1,
        verbose=False
):
    """
    Enriches ASO_df with rank-specific off-target details relative to the specific cell line.

    Instead of summing the Top N genes, this creates columns for the 1st, 2nd, ... Nth
    highest expressed gene individually.

    Returns 3 columns per Rank/Cutoff combo:
      1. Score
      2. Gene Name (String)
      3. Expression Value (Float)
    """
    ASO_df = ASO_df.copy()
    feature_names = []

    # Initialize pandarallel if needed
    if n_cores > 1:
        from pandarallel import pandarallel
        try:
            pandarallel.initialize(nb_workers=n_cores, verbose=0)
        except RuntimeError:
            pass

    # Dictionary to hold lists of partial Series for each new column
    # Structure: accumulator[column_name] = [series_group1, series_group2, ...]
    accumulator = {}

    print(f"Grouping by cell line to calculate specific ranks 1 to {max_rank}...")
    grouped = ASO_df.groupby(CELL_LINE_DEPMAP)

    for cell_line, group_df in grouped:
        if verbose:
            print("Analyzing cell line: ", cell_line)
        # 1. Validation: Check if we have data for this cell line
        if cell_line not in cell_line2data:
            print("Warning! missing cell line: {}".format(cell_line))
            # Missing context: fill this group with NaNs later
            continue

        # 2. Get the specific transcriptome (only up to the max rank needed)
        specific_df = cell_line2data[cell_line].head(max_rank)

        # 3. Identify Normalization Column (Robustness)
        norm_col = next((c for c in specific_df.columns if 'expression_norm' in c), 'expression_norm')

        # Iterate 0 to max_rank-1 (i.e., Rank 1 to Rank N)
        for rank_idx in range(max_rank):
            if verbose:
                print(f"Rank: {rank_idx}")
            current_rank = rank_idx + 1  # 1-based index for naming

            # If the cell line has fewer genes than the requested rank, skip or fill 0
            if rank_idx >= len(specific_df):
                # You might want to log this or handle it.
                # For now, we just won't compute, leaving it as default (0/NaN) initialization.
                continue

            # Get the gene data for this specific rank
            gene_row = specific_df.iloc[rank_idx]
            gene_name = gene_row['Gene'].split()[0]
            if verbose:
                print("Analyzing gene: ", gene_name)
            # Expression value tuple: (TPM, Norm)
            exp_val_tuple = (
                gene_row.get('expression_TPM', 0),
                gene_row.get(norm_col, gene_row.get('expression_norm', 0))
            )
            # We usually save the normalized value in the dataframe column
            exp_val_scalar = exp_val_tuple[1]

            # Validation: Do we have the sequence?
            if gene_name not in gene_to_data:
                continue

            # Construct Single-Gene Maps for the worker
            single_seq_map = {gene_name: gene_to_data[gene_name].full_mrna}
            single_exp_map = {gene_name: exp_val_tuple}

            for cutoff in cutoff_list:
                # define column names
                base_col = f"OT_Spec_Rank{current_rank}_c{cutoff}"
                col_score = f"{base_col}_Score"
                col_gene = f"{base_col}_Gene"
                col_exp = f"{base_col}_Exp"

                # Ensure lists exist in accumulator
                if col_score not in accumulator: accumulator[col_score] = []
                if col_gene not in accumulator:  accumulator[col_gene] = []
                if col_exp not in accumulator:   accumulator[col_exp] = []

                # Compute Score
                # Optimization: Use parallel only for large groups
                if n_cores > 1 and len(group_df) > n_cores:
                    scores = group_df.parallel_apply(
                        compute_single_row,
                        axis=1,
                        args=(single_seq_map, single_exp_map, cutoff, method)
                    )
                else:
                    scores = group_df.apply(
                        compute_single_row,
                        axis=1,
                        args=(single_seq_map, single_exp_map, cutoff, method)
                    )

                # Store Results
                accumulator[col_score].append(scores)

                # Create Series for Metadata (Gene Name and Expression)
                # These are constant for the whole group, but must align with index
                gene_series = pd.Series([gene_name] * len(group_df), index=group_df.index)
                exp_series = pd.Series([exp_val_scalar] * len(group_df), index=group_df.index)

                accumulator[col_gene].append(gene_series)
                accumulator[col_exp].append(exp_series)

                # Track feature names (only add unique ones)
                if col_score not in feature_names:
                    feature_names.extend([col_score, col_gene, col_exp])

    # 4. Reassemble Dataframe
    print("Reassembling chunks...")
    for col_name, chunks in accumulator.items():
        if chunks:
            full_series = pd.concat(chunks)
            # Reindex ensures that rows not processed (e.g. missing cell line data) get NaN/None
            ASO_df[col_name] = full_series.reindex(ASO_df.index)
        else:
            # If no data was computed (empty inputs?), initialize empty
            ASO_df[col_name] = np.nan

    return ASO_df, feature_names