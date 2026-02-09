import numpy as np
from joblib import Parallel, delayed
from scipy.stats import entropy
from tqdm import tqdm

from tauso.data.consts import CELL_LINE_DEPMAP
from tauso.features.rbp.pwm_helper import calculate_total_affinity


def populate_rbp_interaction_features(df, rbp_map, pwm_db, expression_matrix,
                                      sequence_col='flank_sequence_50',
                                      n_jobs=32):
    """
    Calculates the interaction term (Affinity * Expression) for each RBP.

    - Calculates affinity for the sequence.
    - Looks up expression for the cell line.
    - Multiplies them.
    - Assigns directly to df by position (ignoring index).
    """

    flank_param = sequence_col.split('_')[-1]

    df = df.loc[:, ~df.columns.duplicated()]

    sequences = df[sequence_col].fillna("").astype(str).tolist()
    cell_ids = df[CELL_LINE_DEPMAP].tolist()
    n_rows = len(sequences)

    # --- 3. FILTER & PREPARE RBP METADATA ---
    # We only process RBPs that have BOTH a Matrix and (optionally) Expression data.
    # If expression data is missing for an RBP, the result is 0 anyway, so we skip it.

    target_tasks = []

    # Get set of genes present in expression matrix for fast checking
    valid_expr_genes = set(expression_matrix.columns)

    for rbp, matrix_ids in rbp_map.items():
        # 1. Must have a valid PWM
        valid_mids = [m for m in matrix_ids if m in pwm_db]
        if not valid_mids:
            continue

        # 2. Must be in expression matrix (otherwise score is always 0)
        # This is a huge optimization: don't calc affinity if we know expr is 0/missing.
        if rbp not in valid_expr_genes:
            continue

        target_tasks.append({
            'name': rbp,
            'matrix': pwm_db[valid_mids[0]],  # Use the first valid matrix
            'expr_gene': rbp,
            'col_name': f"RBP_{rbp}_expr_aff_{flank_param}"
        })

    if not target_tasks:
        print("Warning: No valid RBP tasks (intersection of PWMs and Expression Matrix).")
        return df, []

    print(f"Calculating interaction features for {len(target_tasks)} RBPs on {n_rows} rows...")
    print(f"Optimization: Calculating affinity only when Expression > 0")

    # --- 4. PARALLEL WORKER ---
    # We process chunks of ROWS.
    # Inside each chunk, we calculate features for ALL RBPs.

    def process_chunk(chunk_seqs, chunk_cells):
        # Result dict: { 'RBP_X_...': [score, score, ...] }
        chunk_results = {task['col_name']: [] for task in target_tasks}

        # We pre-fetch the expression vectors for this chunk of cells to avoid repeated DF lookups
        # subset_expr: DataFrame of shape (len(chunk), len(target_tasks))
        # reindex handles missing cell IDs by putting NaNs (which we fill with 0)
        subset_expr = expression_matrix.reindex(chunk_cells, fill_value=0.0)

        # Iterate over tasks (RBPs)
        for task in target_tasks:
            rbp_name = task['name']
            matrix = task['matrix']
            col_name = task['col_name']

            # Get the expression vector for this RBP for this chunk
            # .values gives us a numpy array, faster iteration
            expr_values = subset_expr[task['expr_gene']].fillna(0.0).values

            scores = []
            for seq, expr_val in zip(chunk_seqs, expr_values):
                # OPTIMIZATION: If expression is 0, interaction is 0.
                # Skip the heavy affinity calculation.
                if expr_val == 0:
                    scores.append(0.0)
                else:
                    aff = calculate_total_affinity(seq, matrix)
                    scores.append(aff * expr_val)

            chunk_results[col_name] = scores

        return chunk_results

    # --- 5. EXECUTION ---
    # Chunking
    chunk_size = max(1, n_rows // (n_jobs * 2))
    chunk_indices = list(range(0, n_rows, chunk_size))

    # Prepare batches
    batches = []
    for i in chunk_indices:
        b_seqs = sequences[i: i + chunk_size]
        b_cells = cell_ids[i: i + chunk_size]
        batches.append((b_seqs, b_cells))

    # Run Parallel
    results = Parallel(n_jobs=n_jobs)(
        delayed(process_chunk)(b_seqs, b_cells) for b_seqs, b_cells in tqdm(batches, desc="Computing Features")
    )

    # --- 6. AGGREGATION & ASSIGNMENT (In-Place) ---
    print("Assigning columns...")

    # Initialize lists in a dict for the final columns
    final_columns = {task['col_name']: [] for task in target_tasks}

    # Merge chunk results
    for res in results:
        for key in final_columns:
            final_columns[key].extend(res[key])

    # Assign to DataFrame
    # Since 'final_columns' are simple lists of exact length matching df,
    # pandas assigns them by position, ignoring the index.
    new_col_names = []
    for col_name, values in final_columns.items():
        df[col_name] = values
        new_col_names.append(col_name)

    print(f"Done. Added {len(new_col_names)} interaction features.")
    return df, new_col_names


def populate_functional_features(df, rbp_map, flank_size=50):
    """
    Aggregates individual RBP interaction columns into functional roles
    and calculates the 'Tug-of-War' (Contrast) features.

    Args:
        df: DataFrame containing individual RBP columns (e.g., RBP_ELAVL1_expr_aff_50)
        rbp_map: Dictionary mapping Gene Name -> Role
        flank_size: The window size suffix (e.g., "50")

    Returns:
        df: Updated DataFrame with new columns
        new_features: List of new column names
    """
    new_features = []

    # Store the column names for the two roles so we can contrast them later
    role_cols = {'stabilizer': None, 'destabilizer': None}

    # --- Step A: Calculate Aggregated Sums ---
    for role in ['stabilizer', 'destabilizer']:
        # 1. Identify which genes in the map belong to this role
        target_genes = [gene for gene, r in rbp_map.items() if r == role]

        # 2. Find the corresponding columns in the DataFrame
        #    (We check if they exist, in case some genes were missing from expression data)
        relevant_columns = [
            f"RBP_{gene}_expr_aff_{flank_size}"
            for gene in target_genes
            if f"RBP_{gene}_expr_aff_{flank_size}" in df.columns
        ]

        if not relevant_columns:
            print(f"Warning: No valid columns found for strict role '{role}'")
            continue

        # 3. Create the Feature: Sum of Interactions
        feature_name = f"RBP_interaction_{role}_{flank_size}"
        df[feature_name] = df[relevant_columns].sum(axis=1)

        new_features.append(feature_name)
        role_cols[role] = feature_name

    # --- Step B: Calculate Contrasts (The "Tug-of-War") ---
    col_stab = role_cols['stabilizer']
    col_destab = role_cols['destabilizer']

    if col_stab and col_destab:
        # 1. Delta: Net Stability (Positive = Stable, Negative = Unstable)
        delta_name = f"RBP_interaction_delta_stab_vs_destab_{flank_size}"
        df[delta_name] = df[col_stab] - df[col_destab]
        new_features.append(delta_name)

        # 2. Ratio: Relative Stability
        #    (We add 1e-6 to denominator to avoid division by zero)
        ratio_name = f"RBP_interaction_ratio_stab_vs_destab_{flank_size}"
        df[ratio_name] = df[col_stab] / (df[col_destab] + 1e-6)
        new_features.append(ratio_name)

    return df, new_features

def populate_complexity_features(df, feature_cols, suffix):
    """
    Calculates role-independent global features:
    1. Total Interaction Load (Sum)
    2. Global Diversity (Entropy)
    """
    # Create a matrix of the relevant columns (add epsilon to avoid div/0)
    # Filling NaNs with 0 is crucial if some lookups failed
    matrix = df[feature_cols].fillna(0.0).values

    # 1. Total Interaction (Sum of all RBPs)
    total_col = f"RBP_interaction_total_{suffix}"
    total_scores = matrix.sum(axis=1)
    df[total_col] = total_scores

    # 2. Global Diversity (Shannon Entropy)
    # Normalize rows to sum to 1 to treat as probabilities
    # Avoid division by zero for rows with 0 interaction
    with np.errstate(divide='ignore', invalid='ignore'):
        probs = matrix / total_scores[:, None]
        # Replace NaNs (from 0/0) with 0
        probs = np.nan_to_num(probs)

    div_col = f"RBP_diversity_global_{suffix}"
    # Calculate entropy (base e by default)
    df[div_col] = entropy(probs, axis=1)

    return df, [total_col, div_col]

