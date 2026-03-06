import pandas as pd
from numba import njit
from scipy.stats import entropy

import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm

from tauso.data.consts import CELL_LINE_DEPMAP, CANONICAL_GENE
from tauso.features.rbp.RBP_features import add_global_complexity_features, add_strict_functional_features, \
    get_background_probs
from tauso.features.rbp.pwm_helper import calculate_total_affinity



@njit(fastmath=True)
def _calculate_affinity_numba_core(seq_indices, pwm_matrix, background_probs):
    seq_len = len(seq_indices)
    motif_len = pwm_matrix.shape[0]

    if seq_len < motif_len:
        return 0.0

    # Numba vectorizes this numpy math perfectly
    weights = np.log2((pwm_matrix + 1e-9) / background_probs)

    total_score = 0.0

    # The classic double loop - Numba excels at this, compiling it down to C
    for i in range(seq_len - motif_len + 1):
        score = 0.0
        valid_window = True

        for pos in range(motif_len):
            base_idx = seq_indices[i + pos]
            if base_idx == -1:
                valid_window = False
                break  # Instantly skip to the next window

            score += weights[pos, base_idx]

        if valid_window and score > 0:
            total_score += score

    return total_score

def calculate_total_affinity_numba(sequence, pwm_matrix, background_probs=None, debug=False):
    # 1. Type guard the matrix
    if hasattr(pwm_matrix, 'values'):
        pwm_matrix = pwm_matrix.values
    else:
        pwm_matrix = np.asarray(pwm_matrix)

    # 2. Handle empty/NaN sequences
    if pd.isna(sequence) or sequence is None:
        return 0.0

    seq_str = str(sequence)
    if len(seq_str) == 0:
        return 0.0

    # 3. Handle background probs
    if background_probs is None:
        background_probs = np.array([0.25, 0.25, 0.25, 0.25], dtype=np.float64)
    else:
        background_probs = np.asarray(background_probs, dtype=np.float64)

    # 4. Map string to integers (Pure Python is fine here before handing off to Numba)
    base_map = {'A': 0, 'C': 1, 'G': 2, 'U': 3, 'T': 3}
    seq_indices = np.array([base_map.get(base, -1) for base in seq_str.upper()], dtype=np.int8)

    # 5. Execute the compiled loop
    # We enforce float64 to ensure math precision matches your original numpy logic perfectly
    return _calculate_affinity_numba_core(seq_indices, pwm_matrix.astype(np.float64), background_probs)



def process_rbp(task, sequences, background_probs_arr):
    """
    Worker function: Processes ONE RBP for ALL sequences.
    """
    matrix = task['matrix']
    col_name = task['col_name']
    n_rows = len(sequences)

    # 1. Pre-allocate a NumPy array (way faster than list.append)
    scores = np.zeros(n_rows, dtype=np.float32)

    # 2. Iterate through sequences
    for i in range(n_rows):
        scores[i] = calculate_total_affinity_numba(sequences[i], matrix, background_probs_arr[i])

    return col_name, scores


def populate_rbp_affinity_features(df, rbp_map, pwm_db, gene_to_data,
                                   sequence_col='flank_sequence_50',
                                   n_jobs=32):
    """
    Calculates the raw Affinity for each RBP (regardless of Expression).
    """
    flank_param = sequence_col.split('_')[-1]
    df = df.loc[:, ~df.columns.duplicated()].copy()  # Use .copy() to avoid SettingWithCopy warnings later

    # Keep sequences as a list for now, as strings don't map well in NumPy
    sequences = df[sequence_col].fillna("").astype(str).tolist()
    n_rows = len(sequences)

    # --- 1. OPTIMIZATION: Create a 2D NumPy array for Background Probs ---
    # This allows Joblib to share memory across workers instead of pickling 180k tiny lists
    default_bg = np.array([0.25, 0.25, 0.25, 0.25], dtype=np.float32)
    bg_list = [
        get_background_probs(gene_to_data[g].full_mrna) if g in gene_to_data else default_bg
        for g in df[CANONICAL_GENE]
    ]
    background_probs_arr = np.array(bg_list, dtype=np.float32)

    # --- 2. FILTER & PREPARE RBP METADATA ---
    target_tasks = []
    for rbp, matrix_ids in rbp_map.items():
        valid_mids = [m for m in matrix_ids if m in pwm_db]
        if not valid_mids:
            continue

        target_tasks.append({
            'name': rbp,
            'matrix': pwm_db[valid_mids[0]],
            'col_name': f"RBP_{rbp}_aff_{flank_param}"
        })

    if not target_tasks:
        print("Warning: No valid RBP tasks found in PWM DB.")
        return df, []

    print(f"Calculating affinity features for {len(target_tasks)} RBPs on {n_rows} rows...")

    # --- 3. EXECUTION: Parallelize over RBPs, not Rows ---
    # joblib uses 'loky' backend by default, which excels at memory mapping large arrays (background_probs_arr)
    results = Parallel(n_jobs=n_jobs)(
        delayed(process_rbp)(task, sequences, background_probs_arr)
        for task in tqdm(target_tasks, desc="Computing RBPs")
    )

    # --- 4. AGGREGATION & ASSIGNMENT ---
    print("Assigning columns...")
    new_col_names = []

    # Because we computed entire columns at a time, assignment is instantaneous
    for col_name, scores in results:
        df[col_name] = scores
        new_col_names.append(col_name)

    print(f"Done. Added {len(new_col_names)} affinity features.")
    return df, new_col_names


def populate_rbp_interaction_features(df, rbp_map, pwm_db, expression_matrix, gene_to_data,
                                      sequence_col='flank_sequence_50',
                                      n_jobs=32):
    """
    Calculates RBP interaction features (Affinity * Expression).
    Includes loud failure checks for common data mismatch issues.
    """
    flank_param = sequence_col.split('_')[-1]

    # ==========================================
    # 1. CRITICAL DATA INTEGRITY CHECKS (LOUD FAILURES)
    # ==========================================

    # CHECK A: Empty Sequences
    # If the sequence column is all empty strings or NaNs, affinity will always be 0.
    if sequence_col not in df.columns:
        raise ValueError(f"❌ MISSING COLUMN: '{sequence_col}' not found in DataFrame.")

    non_empty_seqs = df[sequence_col].astype(str).str.len() > 0
    if not non_empty_seqs.any():
        raise ValueError(
            f"❌ DATA ERROR: Column '{sequence_col}' contains only empty strings or NaNs.\n"
            "   Fix your 'add_external_mrna_and_context_columns' step first!"
        )

    # CHECK B: Row Mismatch (Cell IDs)
    # Force string types to ensure "123" matches "123"
    df[CELL_LINE_DEPMAP] = df[CELL_LINE_DEPMAP].astype(str)
    expression_matrix.index = expression_matrix.index.astype(str)

    common_cells = set(df[CELL_LINE_DEPMAP]).intersection(expression_matrix.index)
    if len(common_cells) == 0:
        raise ValueError(
            "❌ DATA ERROR: No overlap between DataFrame Cell IDs and Expression Matrix index.\n"
            f"   DF Sample: {df[CELL_LINE_DEPMAP].iloc[:3].tolist()}\n"
            f"   Matrix Sample: {expression_matrix.index[:3].tolist()}"
        )

    # CHECK C: Column Mismatch (Gene Names)
    # Check if any RBP from the map exists in the expression matrix
    valid_expr_genes = set(expression_matrix.columns)
    valid_rbps = [r for r in rbp_map.keys() if r in valid_expr_genes]

    if not valid_rbps:
        raise ValueError(
            "❌ DATA ERROR: None of the RBPs in 'rbp_map' exist in the Expression Matrix columns.\n"
            "   Check gene name formats (e.g. 'MYC' vs 'MYC (4609)')."
        )

    print(f"{len(common_cells)} cells, {len(valid_rbps)} valid RBPs found.")
    # ==========================================

    # --- 2. PREPARATION ---
    df = df.loc[:, ~df.columns.duplicated()]

    # Convert to standard types for pickle/joblib safety
    sequences = df[sequence_col].fillna("").astype(str).tolist()
    cell_ids = df[CELL_LINE_DEPMAP].tolist()

    # Pre-calculate background probs
    default_bg = np.array([0.25, 0.25, 0.25, 0.25])
    background_probs = []
    for g in df[CANONICAL_GENE]:
        if g in gene_to_data:
            background_probs.append(get_background_probs(gene_to_data[g].full_mrna))
        else:
            print("Using default probabilities for gene: ", g)
            background_probs.append(default_bg)

    n_rows = len(sequences)

    # --- 3. FILTER TASKS ---
    target_tasks = []
    for rbp, matrix_ids in rbp_map.items():
        # Must be in expression matrix (already checked generally, now specific filter)
        if rbp not in valid_expr_genes:
            continue

        # Must have valid PWM
        valid_mids = [m for m in matrix_ids if m in pwm_db]
        if not valid_mids:
            continue

        target_tasks.append({
            'name': rbp,
            'matrix': pwm_db[valid_mids[0]],
            'expr_gene': rbp,
            'col_name': f"RBP_{rbp}_expr_aff_{flank_param}"
        })

    if not target_tasks:
        print("Warning: No valid RBP tasks remaining after filtering.")
        return df, []

    print(f"Calculating interaction features for {len(target_tasks)} RBPs on {n_rows} rows...")

    # --- 4. PARALLEL WORKER ---
    def process_chunk(chunk_seqs, chunk_cells, chunk_background_probs):
        chunk_results = {task['col_name']: [] for task in target_tasks}

        # Optimized reindex: gets all necessary expression data for this chunk at once
        subset_expr = expression_matrix.reindex(chunk_cells, fill_value=np.nan)

        for task in target_tasks:
            matrix = task['matrix']
            col_name = task['col_name']

            # Extract expression column (fast numpy access)
            expr_values = subset_expr[task['expr_gene']].values

            scores = []
            for seq, expr_val, bg_probs in zip(chunk_seqs, expr_values, chunk_background_probs):
                # OPTIMIZATION: Handle missing data and true zeros
                if pd.isna(expr_val):
                    scores.append(np.nan)
                elif expr_val == 0:
                    scores.append(0.0)
                else:
                    # Only pay the cost of affinity calculation if expression > 0
                    aff = calculate_total_affinity_numba(seq, matrix, bg_probs)
                    scores.append(aff * expr_val)

            chunk_results[col_name] = scores

        return chunk_results

    # --- 5. EXECUTION ---
    chunk_size = max(1, n_rows // (n_jobs * 2))
    chunk_indices = list(range(0, n_rows, chunk_size))

    batches = []
    for i in chunk_indices:
        b_seqs = sequences[i: i + chunk_size]
        b_cells = cell_ids[i: i + chunk_size]
        b_background_probs = background_probs[i: i + chunk_size]
        batches.append((b_seqs, b_cells, b_background_probs))

    results = Parallel(n_jobs=n_jobs)(
        delayed(process_chunk)(b_seqs, b_cells, b_bg)
        for b_seqs, b_cells, b_bg in tqdm(batches, desc="Computing Batches")
    )

    # --- 6. AGGREGATION ---
    print("Assigning columns...")
    final_columns = {task['col_name']: [] for task in target_tasks}

    for res in results:
        for key in final_columns:
            final_columns[key].extend(res[key])

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

def populate_complexity_features(df, feature_cols, suffix, type='generic'):
    """
    Calculates role-independent global features:
    1. Total Interaction Load (Sum)
    2. Global Diversity (Entropy)
    """
    # Create a matrix of the relevant columns (add epsilon to avoid div/0)
    # Filling NaNs with 0 is crucial if some lookups failed
    matrix = df[feature_cols].fillna(0.0).values

    # 1. Total Interaction (Sum of all RBPs)
    total_col = f"RBP_interaction_total_{suffix}_{type}"
    total_scores = matrix.sum(axis=1)
    df[total_col] = total_scores

    # 2. Global Diversity (Shannon Entropy)
    # Normalize rows to sum to 1 to treat as probabilities
    # Avoid division by zero for rows with 0 interaction
    with np.errstate(divide='ignore', invalid='ignore'):
        probs = matrix / total_scores[:, None]
        # Replace NaNs (from 0/0) with 0
        probs = np.nan_to_num(probs)

    div_col = f"RBP_diversity_global_{suffix}_{type}"
    # Calculate entropy (base e by default)
    df[div_col] = entropy(probs, axis=1)

    return df, [total_col, div_col]


def populate_rbp_region(df, region_name, rbp_map, pwm_db, expr_matrix, role_map, gene_to_data, n_jobs=32):
    """
    1. Generates raw RBP features for a specific region (Left/Core/Right).
    2. Calculates aggregate Complexity & Functional scores.
    3. Drops the raw columns to return a clean DataFrame.
    """
    print(f"\n--- Processing Region: {region_name.upper()} ---")

    # A. Run the Engine (Generate Raw Features)
    # We assume columns like 'seq_region_left' already exist
    df, raw_feats = populate_rbp_interaction_features(
        df, rbp_map, pwm_db, expr_matrix, gene_to_data,
        sequence_col=f"seq_region_{region_name}",
        n_jobs=n_jobs
    )

    # B. Calculate Aggregates
    df, complexity_feats = add_global_complexity_features(df, raw_feats, suffix=region_name)
    df, functional_feats = add_strict_functional_features(df, role_map, suffix=region_name)

    # C. Cleanup (Drop raw columns immediately)
    df.drop(columns=raw_feats, inplace=True)
    print(f"Dropped {len(raw_feats)} raw columns for {region_name}.")

    # Return DF and list of new columns to save
    return df, complexity_feats + functional_feats
