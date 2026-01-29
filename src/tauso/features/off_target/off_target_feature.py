from .add_off_target_feat import  compute_general_vector_and_cache, compute_specific_vector_optimized
import pandas as pd


def add_depmap_id_column(ASO_df, name_to_depmap):
    ASO_df = ASO_df.copy()
    ASO_df["DepMap_ID"] = ASO_df["Cell_line"].map(name_to_depmap)
    return ASO_df

def compute_off_target_vector(
        ASO_df,
        gene_to_data,
        cell_line2data,
        cutoff,
        method,
        top_n):
    """
    Computes two vectors:
    1. specific_scores: {index: score} based on the specific cell line's transcriptome.
    2. general_scores:  {index: score} based on the 'general' transcriptome.

    Optimization: Caches energies from the General calculation to speed up Specific calculations.
    """

    # 1. Initialize Outputs
    specific_scores = {}
    general_scores = {}

    print("Phase 1: Computing General Scores and Caching Energies...")
    # 3. Retrieve General Transcriptome
    if 'general' not in cell_line2data:
        raise ValueError("Key 'general' not found in cell_line2data dictionary.")

    # Assuming the structure is {key: [df]}, extract the dataframe
    general_df = cell_line2data['general'].head(top_n)

    # 4. Phase 1: Compute General Baseline (Returns scores AND the energy cache)
    # energy_cache structure: { ASO_Sequence: { Gene_Name: (Energy, Transcript_Sequence) } }
    general_scores, energy_cache = compute_general_vector_and_cache(
        ASO_df=ASO_df,
        gene_to_data=gene_to_data,
        general_df=general_df,
        cutoff=cutoff,
        method=method
    )

    print("Phase 2: Computing Specific Scores using Baseline + Delta...")

    # 5. Phase 2: Compute Specific Scores (Grouped by Cell Line)
    for depmap_id, group_df in ASO_df.groupby("DepMap_ID"):

        if pd.isna(depmap_id) or depmap_id not in cell_line2data:
            continue

        # Get Specific Transcriptome
        transcriptome_df = cell_line2data[depmap_id].head(top_n)

        # Calculate Specific scores using the Energy Cache
        result = compute_specific_vector_optimized(
            specific_df=transcriptome_df,
            gene_to_data=gene_to_data,
            ASO_df=group_df,
            energy_cache=energy_cache,
            cutoff=cutoff,
            method=method
        )
        specific_scores.update(result)

    return specific_scores, general_scores


def serialize_feature_name(method, top_n, cutoff, is_specific):
    if is_specific:
        return f"off_target_score_specific_{method}_n{top_n}_c{cutoff}"
    return f"off_target_score_general_{method}_n{top_n}_c{cutoff}"

def run_off_target_pipeline(
        ASO_df,
        gene_to_data,
        cell_line2data,
        top_n_list,
        cutoff_list,
        method
):
    """
    Enriches ASO_df with off-target scores instead of saving to files.

    Returns:
        tuple: (enriched_ASO_df, feature_names_list)
    """
    # Create a copy to ensure we don't modify the original dataframe slice in place
    # or trigger SettingWithCopy warnings.
    ASO_df = ASO_df.copy()

    feature_names = []

    for top_n in top_n_list:
        for cutoff in cutoff_list:
            print(f"Running config: TopN={top_n}, Cutoff={cutoff}")

            # Call the computation function
            spec_scores, gen_scores = compute_off_target_vector(
                ASO_df=ASO_df,
                gene_to_data=gene_to_data,
                cell_line2data=cell_line2data,
                cutoff=cutoff,
                method=method,
                top_n=top_n
            )

            # --- Process Specific Scores ---
            spec_col_name = serialize_feature_name(method, top_n, cutoff, is_specific=True)

            # Map the scores dictionary to the dataframe index
            # This assumes spec_scores keys correspond to ASO_df indices
            ASO_df[spec_col_name] = ASO_df.index.map(spec_scores)
            feature_names.append(spec_col_name)

            # --- Process General Scores ---
            if gen_scores:
                gen_col_name = serialize_feature_name(method, top_n, cutoff, is_specific=False)

                # Map the scores dictionary to the dataframe index
                ASO_df[gen_col_name] = ASO_df.index.map(gen_scores)
                feature_names.append(gen_col_name)

    return ASO_df, feature_names


def run_off_target_pipeline_general(
        ASO_df,
        gene_to_data,
        cell_line2data,
        top_n_list,
        cutoff_list,
        method
):
    """
    Enriches ASO_df with off-target scores instead of saving to files.

    Returns:
        tuple: (enriched_ASO_df, feature_names_list)
    """
    # Create a copy to ensure we don't modify the original dataframe slice in place
    # or trigger SettingWithCopy warnings.
    ASO_df = ASO_df.copy()

    feature_names = []

    if 'general' not in cell_line2data:
        raise ValueError("Key 'general' not found in cell_line2data dictionary.")
    general_df_all = cell_line2data['general']

    for top_n in top_n_list:
        # Assuming the structure is {key: [df]}, extract the dataframe
        general_df = general_df_all.head(top_n)

        for cutoff in cutoff_list:
            print(f"Running config: TopN={top_n}, Cutoff={cutoff}")
            # Call the computation function
            gen_scores, _ = compute_general_vector_and_cache(
                ASO_df=ASO_df,
                gene_to_data=gene_to_data,
                general_df=general_df,
                cutoff=cutoff,
                method=method
            )

            gen_col_name = serialize_feature_name(method, top_n, cutoff, is_specific=False)

            # Map the scores dictionary to the dataframe index
            ASO_df[gen_col_name] = ASO_df.index.map(gen_scores)
            feature_names.append(gen_col_name)

    return ASO_df, feature_names