import pandas as pd

from ...common.gtf import filter_gtf_genes

def general_exp_data(EXP_path, db, filter_mode="protein_coding"):
    """
    Loads expression data, filters for valid genes based on GTF (e.g., protein_coding),
    averages across cell lines, and converts to TPM.
    """
    # 1. Get the set of allowed gene names (the "clean" list)
    allowed_names = filter_gtf_genes(db, filter_mode)

    # 2. Load expression matrix
    exp_data = pd.read_csv(EXP_path)

    # 3. Filter the columns
    # First, identify columns that look like genes: "GeneName (ID)"
    potential_gene_cols = [c for c in exp_data.columns if "(" in c]

    valid_cols = []
    for col in potential_gene_cols:
        # Extract the gene name part: "TP53 (7157)" -> "TP53"
        # We split by " (" and take the first part
        gene_name_extracted = col.split(" (")[0]

        # Check if this gene is in our allowed whitelist
        if gene_name_extracted in allowed_names:
            valid_cols.append(col)

    # Apply the filter to the dataframe
    print(f"Filtering: Kept {len(valid_cols)} genes out of {len(potential_gene_cols)} total columns.")
    exp_data = exp_data[valid_cols]

    # 4. Average across cell lines (axis=0) → one value per gene
    mean_exp = exp_data.mean(axis=0)

    # 5. Format the output
    mean_exp_data = mean_exp.reset_index()
    mean_exp_data.columns = ["Gene", "general_expression_norm"]

    # 6. Back-transform to TPM (assuming input was log2(TPM+1) or similar)
    mean_exp_data["expression_TPM"] = 2 ** mean_exp_data["general_expression_norm"]

    # Optional: Clean up the "Gene" column in the output to remove the "(ID)"
    # so it matches your other data formats if needed.
    # mean_exp_data["Gene"] = mean_exp_data["Gene"].apply(lambda x: x.split(" (")[0])

    return mean_exp_data.sort_values("general_expression_norm", ascending=False)
