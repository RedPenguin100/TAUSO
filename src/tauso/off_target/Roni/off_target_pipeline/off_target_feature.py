from src.tauso.off_target.Roni.off_target_pipeline.add_off_target_feat import  add_off_target_feat
import pandas as pd
import os


def add_depmap_id_column(ASO_df, name_to_depmap):
    ASO_df = ASO_df.copy()
    ASO_df["DepMap_ID"] = ASO_df["Cell_line"].map(name_to_depmap)
    return ASO_df

def compute_off_target_vector(
    ASO_df,
    cell_line2data,
    cutoff,
    method,
    top_n,
    name_to_depmap
):
    """
    Compute off-target scores for a set of ASOs.

    Returns
    -------
    dict
        Mapping ASO index -> off-target score
    """
    scores = {}

    ASO_df = add_depmap_id_column(ASO_df, name_to_depmap)

    for depmap_id, group_df in ASO_df.groupby("DepMap_ID"):

        if pd.isna(depmap_id):
            continue

        if depmap_id not in cell_line2data:
            continue

        transcriptome_df = cell_line2data[depmap_id][-1].head(top_n)
        cell_line2df = {depmap_id: transcriptome_df}

        result = add_off_target_feat(
            cell_line_dict=name_to_depmap,
            cell_line2df=cell_line2df,
            ASO_df=group_df,
            cutoff=cutoff,
            method=method
        )
        scores.update(result)

    return scores

def run_single_configuration(
    ASO_df,
    cell_line2data,
    output_dir,
    top_n,
    cutoff,
    method,
    method_label
):
    scores = compute_off_target_vector(
        ASO_df=ASO_df,
        cell_line2data=cell_line2data,
        cutoff=cutoff,
        method=method,
        top_n=top_n
    )

    df = pd.DataFrame({
        "index": list(scores.keys()),
        f"off_target_score_{method_label}": list(scores.values())
    })

    out_name = f"off_target.{method_label}.top{top_n}.cutoff{cutoff}.csv"
    df.to_csv(os.path.join(output_dir, out_name), index=False)


def run_off_target_pipeline(
    ASO_df,
    cell_line2data,
    output_dir,
    top_n_list,
    cutoff_list,
    method,
    name_to_depmap
):
    """
    Run off-target analysis across multiple configurations.
    """
    method_dict = {
        0: "ARTM_log",
        1: "GEO",
        2: "ARTM",
        3: "MS",
        4: "MECH"
    }
    method_label = method_dict[method]

    os.makedirs(output_dir, exist_ok=True)

    for top_n in top_n_list:
        for cutoff in cutoff_list:
            scores = compute_off_target_vector(
                ASO_df=ASO_df,
                cell_line2data=cell_line2data,
                cutoff=cutoff,
                method=method,
                top_n=top_n,
                name_to_depmap=name_to_depmap
            )

            df = pd.DataFrame({
                "index": list(scores.keys()),
                f"off_target_score_{method_label}": list(scores.values())
            })

            out_name = f"off_target.{method_label}.top{top_n}.cutoff{cutoff}.csv"
            df.to_csv(os.path.join(output_dir, out_name), index=False)
