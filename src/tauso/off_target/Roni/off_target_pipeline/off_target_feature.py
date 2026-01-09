# from src.tauso.off_target.Roni.off_target_pipeline.add_off_target_feat import  add_off_target_feat
from src.tauso.off_target.Roni.off_target_pipeline.add_off_target_feat import  compute_general_vector_and_cache, compute_specific_vector_optimized
import pandas as pd
import os


def add_depmap_id_column(ASO_df, name_to_depmap):
    ASO_df = ASO_df.copy()
    ASO_df["DepMap_ID"] = ASO_df["Cell_line"].map(name_to_depmap)
    return ASO_df

# def compute_off_target_vector(
#     ASO_df,
#     cell_line2data,
#     cutoff,
#     method,
#     top_n,
#     name_to_depmap
# ):
#     """
#     Compute off-target scores for a set of ASOs.
#
#     Returns
#     -------
#     dict
#         Mapping ASO index -> off-target score
#     """
#     scores = {}
#
#     ASO_df = add_depmap_id_column(ASO_df, name_to_depmap)
#
#     for depmap_id, group_df in ASO_df.groupby("DepMap_ID"):
#
#         if pd.isna(depmap_id):
#             continue
#
#         if depmap_id not in cell_line2data:
#             continue
#
#         transcriptome_df = cell_line2data[depmap_id][-1].head(top_n)
#         cell_line2df = {depmap_id: transcriptome_df}
#
#         result = add_off_target_feat(
#             cell_line_dict=name_to_depmap,
#             cell_line2df=cell_line2df,
#             ASO_df=group_df,
#             cutoff=cutoff,
#             method=method
#         )
#         scores.update(result)
#
#     return scores


def compute_off_target_vector(
        ASO_df,
        cell_line2data,
        cutoff,
        method,
        top_n,
        name_to_depmap
):
    """
    Computes two vectors:
    1. specific_scores: {index: score} based on the specific cell line's transcriptome.
    2. general_scores:  {index: score} based on the 'general' transcriptome.

    Optimization: Caches energies from the General calculation to speed up Specific calculations.
    """

    # 1. Initialize Outputs
    specific_scores = {}
    general_scores = {}

    # 2. Add DepMap IDs
    ASO_df = add_depmap_id_column(ASO_df, name_to_depmap)

    # 3. Retrieve General Transcriptome
    if 'general' not in cell_line2data:
        raise ValueError("Key 'general' not found in cell_line2data dictionary.")

    # Assuming the structure is {key: [df]}, extract the dataframe
    general_df = cell_line2data['general'][0].head(top_n)

    print("Phase 1: Computing General Scores and Caching Energies...")

    # 4. Phase 1: Compute General Baseline (Returns scores AND the energy cache)
    # energy_cache structure: { ASO_Sequence: { Gene_Name: (Energy, Transcript_Sequence) } }
    general_scores, energy_cache = compute_general_vector_and_cache(
        ASO_df=ASO_df,
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
        transcriptome_df = cell_line2data[depmap_id][-1].head(top_n)

        # Calculate Specific scores using the Energy Cache
        result = compute_specific_vector_optimized(
            specific_df=transcriptome_df,
            ASO_df=group_df,
            energy_cache=energy_cache,
            cutoff=cutoff,
            method=method
        )
        specific_scores.update(result)

    return specific_scores, general_scores

# def run_single_configuration(
#     ASO_df,
#     cell_line2data,
#     output_dir,
#     top_n,
#     cutoff,
#     method,
#     method_label
# ):
#     scores = compute_off_target_vector(
#         ASO_df=ASO_df,
#         cell_line2data=cell_line2data,
#         cutoff=cutoff,
#         method=method,
#         top_n=top_n
#     )
#
#     df = pd.DataFrame({
#         "index": list(scores.keys()),
#         f"off_target_score_{method_label}": list(scores.values())
#     })
#
#     out_name = f"off_target.{method_label}.top{top_n}.cutoff{cutoff}.csv"
#     df.to_csv(os.path.join(output_dir, out_name), index=False)

def run_single_configuration(
        ASO_df,
        cell_line2data,
        output_dir,
        top_n,
        cutoff,
        method,
        method_label,
        name_to_depmap  # <--- Added this argument
):
    # 1. Call the compute function and unpack the two result dictionaries
    spec_scores, gen_scores = compute_off_target_vector(
        ASO_df=ASO_df,
        cell_line2data=cell_line2data,
        cutoff=cutoff,
        method=method,
        top_n=top_n,
        name_to_depmap=name_to_depmap
    )

    # 2. Save Specific Scores (Always expected)
    df_spec = pd.DataFrame({
        "index": list(spec_scores.keys()),
        f"off_target_score_{method_label}_specific": list(spec_scores.values())
    })

    out_name_spec = f"off_target.{method_label}.specific.top{top_n}.cutoff{cutoff}.csv"
    df_spec.to_csv(os.path.join(output_dir, out_name_spec), index=False)
    print(f"Saved Specific scores to {out_name_spec}")

    # 3. Save General Scores (Only if they exist)
    if gen_scores:
        df_gen = pd.DataFrame({
            "index": list(gen_scores.keys()),
            f"off_target_score_{method_label}_general": list(gen_scores.values())
        })

        out_name_gen = f"off_target.{method_label}.general.top{top_n}.cutoff{cutoff}.csv"
        df_gen.to_csv(os.path.join(output_dir, out_name_gen), index=False)
        print(f"Saved General scores to {out_name_gen}")


def run_off_target_pipeline(
        ASO_df,
        cell_line2data,
        output_dir,
        top_n_list,
        cutoff_list,
        method,
        name_to_depmap
):
    method_dict = {0: "ARTM_log", 1: "GEO", 2: "ARTM", 3: "MS", 4: "MECH"}
    method_label = method_dict.get(method, "Unknown")

    os.makedirs(output_dir, exist_ok=True)

    for top_n in top_n_list:
        for cutoff in cutoff_list:
            print(f"Running config: TopN={top_n}, Cutoff={cutoff}")

            # This calls the function above
            spec_scores, gen_scores = compute_off_target_vector(
                ASO_df=ASO_df,
                cell_line2data=cell_line2data,
                cutoff=cutoff,
                method=method,
                top_n=top_n,
                name_to_depmap=name_to_depmap
            )

            # Save Specific Scores
            df_spec = pd.DataFrame({
                "index": list(spec_scores.keys()),
                f"off_target_score_{method_label}_specific": list(spec_scores.values())
            })
            out_spec = f"off_target.{method_label}.specific.top{top_n}.cutoff{cutoff}.csv"
            df_spec.to_csv(os.path.join(output_dir, out_spec), index=False)

            # Save General Scores
            if gen_scores:
                df_gen = pd.DataFrame({
                    "index": list(gen_scores.keys()),
                    f"off_target_score_{method_label}_general": list(gen_scores.values())
                })
                out_gen = f"off_target.{method_label}.general.top{top_n}.cutoff{cutoff}.csv"
                df_gen.to_csv(os.path.join(output_dir, out_gen), index=False)
