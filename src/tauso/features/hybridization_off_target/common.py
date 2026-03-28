import pandas as pd


# TODO: move this function to a folder / file about expression


def get_general_expression_of_genes(EXP_path, valid_genes):
    """
    Loads expression data, filters for valid genes based on GTF (e.g., protein_coding),
    averages across cell lines, and converts to TPM.
    """
    # OPTIMIZATION 1: Convert list to set for O(1) instant lookups
    valid_genes_set = set(valid_genes)

    # OPTIMIZATION 2: Read ONLY the first row (header) to get column names instantly
    header_df = pd.read_csv(EXP_path, nrows=0)

    # Identify columns that look like genes: "GeneName (ID)"
    potential_gene_cols = [c for c in header_df.columns if "(" in c]

    valid_cols = []
    for col in potential_gene_cols:
        # Extract the gene name part: "TP53 (7157)" -> "TP53"
        gene_name_extracted = col.split(" (")[0]

        # Check if this gene is in our allowed whitelist
        if gene_name_extracted in valid_genes_set:
            valid_cols.append(col)

    print(f"Filtering: Kept {len(valid_cols)} genes out of {len(potential_gene_cols)} total columns.")

    # Safety check: if no genes matched, return an empty dataframe with correct columns
    if not valid_cols:
        return pd.DataFrame(columns=["Gene", "expression_norm", "expression_TPM"])

    # OPTIMIZATION 3: Load the massive CSV, but ONLY parse and load the 'valid_cols'
    exp_data = pd.read_csv(EXP_path, usecols=valid_cols)

    # 4. Average across cell lines (axis=0) → one value per gene
    mean_exp = exp_data.mean(axis=0)

    # 5. Format the output
    mean_exp_data = mean_exp.reset_index()
    mean_exp_data.columns = ["Gene", "expression_norm"]

    # 6. Back-transform to TPM (assuming input was log2(TPM+1) or similar)
    mean_exp_data["expression_TPM"] = 2 ** mean_exp_data["expression_norm"]

    # Optional: Clean up the "Gene" column in the output to remove the "(ID)"
    # so it matches your other data formats if needed.
    # mean_exp_data["Gene"] = mean_exp_data["Gene"].apply(lambda x: x.split(" (")[0])

    return mean_exp_data.sort_values("expression_norm", ascending=False)
