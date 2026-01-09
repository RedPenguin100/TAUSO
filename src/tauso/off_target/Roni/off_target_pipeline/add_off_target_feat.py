# from src.tauso.hybridization.fast_hybridization import get_trigger_mfe_scores_by_risearch
from src.tauso.off_target.Roni.external.fold import get_trigger_mfe_scores_by_risearch
from src.tauso.off_target.Roni.off_target_pipeline.off_target_functions import dna_to_rna_reverse_complement, parse_risearch_output, aggregate_off_targets

import pandas as pd
import numpy as np

print("add off target feat")



# path = '/home/oni/ASOdesign/scripts/data_genertion/cell_line_expression/'
# sequences = '.mutated_transcriptome_premRNA.merged.csv'
#
#
#
# # Set base directory relative to script
# script_dir = os.path.dirname(__file__)
# data_dir = os.path.abspath(os.path.join(script_dir, "..", "data_genertion"))
#
# # ========================================== Pre-Processing Main ASO data =============================================
# # Load the main ASO dataset
# data_path = os.path.join(data_dir, "data_updated_inhibition.csv")
# may_df = pd.read_csv(data_path)
# may_df = may_df.head(1)
#
# # Expression files path
#
# # Load each transcriptome file
# A431_df = pd.read_csv(path+celline_list[0]+sequences)
# print("1")
# NCI_H460_df = pd.read_csv(path+celline_list[1]+sequences)
# print("2")
# SH_SY5Y_df = pd.read_csv(path+celline_list[2]+sequences)
# print("3")
# HeLa_df = pd.read_csv(path+celline_list[3]+sequences)
# print("4")
# HepG2_df = pd.read_csv(path+celline_list[4]+sequences)
# print("5")
# U_251MG_df = pd.read_csv(path+celline_list[5]+sequences)
# print("6")
#
# # ============================ Cut to top n expressed ==============================
# n = 1
# A431_df = A431_df.head(n)
# NCI_H460_df = NCI_H460_df.head(n)
# SH_SY5Y_df = SH_SY5Y_df.head(n)
# HeLa_df = HeLa_df.head(n)
# HepG2_df = HepG2_df.head(n)
# U_251MG_df = U_251MG_df.head(n)
#  # =============================================== Cell lines  ========================================================
# cell_line_list = ['ACH-001328', 'ACH-000463', 'ACH-001188', 'ACH-001086', 'ACH-000739', 'ACH-000232']
#
# cell_line2df_ = {'A431':A431_df,
#                 'A-431':A431_df,
#                 'NCI-H460':NCI_H460_df,
#                 'SH_SY5Y':SH_SY5Y_df,
#                 'HeLa':HeLa_df,
#                 'HepG2':HepG2_df,
#                 'U-251MG':U_251MG_df}
#
# =============================================== Function ============================================================

"""
def add_off_target_feat(cell_line_dict, cell_line2df, ASO_df, cutoff, method=0):
    print(cell_line2df.keys())
    index_score_vec = {}
    index_general_vec = {}
    RT = 0.616  # for mechanical statistics

    for _, aso_row in ASO_df.iterrows():

        index = aso_row['index']
        exp_cell_line = aso_row['Cell_line']
        cell_line = cell_line_dict.get(exp_cell_line)
        ASO_seq = aso_row['Sequence']
        target_gene = aso_row['Canonical Gene Name']

        if cell_line not in cell_line2df:
            continue

        trigger = dna_to_rna_reverse_complement(ASO_seq)
        curr_df = cell_line2df[cell_line]

        # ---------------- Collect sequences + expression ----------------
        name_to_seq = {}
        name_to_exp_TPM = {}
        name_to_exp_norm = {}

        for _, gene_row in curr_df.iterrows():
            gene = gene_row['Gene'].split()[0]
            if gene == target_gene:
                continue

            seq = gene_row.get('Mutated Transcript Sequence')
            if pd.isna(seq):
                seq = gene_row.get('Original Transcript Sequence')

            if pd.isna(seq):
                continue

            name_to_seq[gene] = seq
            name_to_exp_TPM[gene] = gene_row['expression_TPM']
            name_to_exp_norm[gene] = gene_row[f'{cell_line}_expression_norm']

        if not name_to_seq:
            continue

        # ---------------- RIsearch ----------------
        result_dict = get_trigger_mfe_scores_by_risearch(
            trigger,
            name_to_seq,
            minimum_score=cutoff,
            parsing_type='2'
        )

        result_df = parse_risearch_output(result_dict)
        result_df_agg = aggregate_off_targets(result_df)

        if result_df_agg.empty:
            index_score_vec[index] = 0.0
            continue

        # ---------------- Scoring ----------------
        score = 0.0

        if method == 0:  # normalized arithmetic mean
            for _, r in result_df_agg.iterrows():
                score += r['energy'] * name_to_exp_norm.get(r['target'], 0)

        elif method == 1:  # TPM arithmetic mean
            for _, r in result_df_agg.iterrows():
                score += r['energy'] * name_to_exp_TPM.get(r['target'], 0)

        elif method == 2:  # weighted arithmetic mean (e^2)
            for _, r in result_df_agg.iterrows():
                g = r['energy']
                e = name_to_exp_TPM.get(r['target'], 0) / 1e6
                if g < 0 and e > 0:
                    score += g * e**2

        elif method == 3:  # geometric mean (log space)
            log_sum = 0.0
            for _, r in result_df_agg.iterrows():
                g = r['energy']
                e = name_to_exp_TPM.get(r['target'], 0) / 1e6
                if g < 0 and e > 0:
                    log_sum += e * np.log(-g)
            score = log_sum

        elif method == 4:  # rank-based decreasing arithmetic weighting
            ranked = result_df_agg.assign(
                exp=[name_to_exp_TPM.get(t, 0) for t in result_df_agg['target']]
            ).sort_values("exp", ascending=False)

            for rank, (_, r) in enumerate(ranked.iterrows(), start=1):
                g = r['energy']
                e = r['exp'] / 1e6
                if g < 0 and e > 0:
                    batch = (rank - 1) // 10
                    weight = 10 / (2 ** batch)
                    score += g * e * weight

        elif method == 5:  # mechanical statistics
            for _, r in result_df_agg.iterrows():
                g = r['energy']
                e = name_to_exp_TPM.get(r['target'], 0) / 1e6
                score += e * np.exp(-g / RT)

        index_score_vec[index] = score

    # normalize mechanical statistics
    if method == 5 and index_score_vec:
        max_val = max(index_score_vec.values())
        if max_val > 0:
            index_score_vec = {k: v / max_val for k, v in index_score_vec.items()}

    return index_score_vec
"""


def compute_general_vector_and_cache(ASO_df, general_df, cutoff, method):
    """
    Calculates the General Vector for all rows and caches the RIsearch energies.
    """
    index_scores = {}
    energy_cache = {}  # Key: ASO_Seq, Value: { Gene: (Energy, Sequence) }

    # 1. Prepare General Data Maps
    general_seq_map = {}
    general_exp_map = {}

    norm_col = next((c for c in general_df.columns if 'expression_norm' in c), 'expression_norm')
    use_norm = (method == 0)

    for _, row in general_df.iterrows():
        gene = row['Gene'].split()[0]
        # General usually uses 'Original', but check both
        seq = row.get('Mutated Transcript Sequence')
        if pd.isna(seq): seq = row.get('Original Transcript Sequence')
        if pd.isna(seq): continue

        general_seq_map[gene] = seq

        # Store expression for scoring
        if use_norm:
            # Fallback to 'expression_norm' if dynamic name fails, else 0
            general_exp_map[gene] = row.get(norm_col, row.get('expression_norm', 0))
        else:
            general_exp_map[gene] = row.get('expression_TPM', 0)

    # 2. Iterate Unique ASO Sequences (Avoid re-calculating identical ASOs)
    unique_asos = ASO_df[['Sequence']].drop_duplicates()

    for _, row in unique_asos.iterrows():
        aso_seq = row['Sequence']
        trigger = dna_to_rna_reverse_complement(aso_seq)

        # --- RIsearch (Heavy Lifting) ---
        result_dict = get_trigger_mfe_scores_by_risearch(
            trigger,
            general_seq_map,
            minimum_score=cutoff,
            parsing_type='2'
        )
        result_df = parse_risearch_output(result_dict)
        result_df_agg = aggregate_off_targets(result_df)

        # Cache Energies AND Sequences {Gene: (Energy, Seq)}
        # We need the Seq to detect mutations in Phase 2
        gene_energy_seq = {}
        if not result_df_agg.empty:
            for _, r in result_df_agg.iterrows():
                g = r['target']
                e = r['energy']
                s = general_seq_map.get(g)
                gene_energy_seq[g] = (e, s)

        energy_cache[aso_seq] = gene_energy_seq

        # Calculate General Score
        # We reconstruct a simple energy dict for the helper function
        simple_energies = {g: e for g, (e, s) in gene_energy_seq.items()}
        score = calculate_score_helper(simple_energies, general_exp_map, method)

        # Assign this score to ALL rows with this ASO sequence
        relevant_indices = ASO_df[ASO_df['Sequence'] == aso_seq]['index']
        for idx in relevant_indices:
            index_scores[idx] = score

    return index_scores, energy_cache


def compute_specific_vector_optimized(specific_df, ASO_df, energy_cache, cutoff, method):
    """
    Calculates scores for a specific cell line.
    OPTIMIZATION: Reuses cached energies if the gene sequence is identical to General.
    """
    index_scores = {}

    # 1. Prepare Specific Data
    specific_exp_map = {}
    specific_seq_map = {}

    norm_col = next((c for c in specific_df.columns if 'expression_norm' in c), 'expression_norm')
    use_norm = (method == 0)

    for _, row in specific_df.iterrows():
        gene = row['Gene'].split()[0]

        # Prefer Mutated, fallback to Original
        seq = row.get('Mutated Transcript Sequence')
        if pd.isna(seq): seq = row.get('Original Transcript Sequence')
        if pd.isna(seq): continue

        specific_seq_map[gene] = seq

        if use_norm:
            specific_exp_map[gene] = row.get(norm_col, row.get('expression_norm', 0))
        else:
            specific_exp_map[gene] = row.get('expression_TPM', 0)

    # 2. Iterate ASOs
    for _, row in ASO_df.iterrows():
        index = row['index']
        aso_seq = row['Sequence']
        target_gene = row['Canonical Gene Name']

        # Retrieve Cache: {Gene: (Energy, General_Seq)}
        cached_data = energy_cache.get(aso_seq, {})

        # 3. Identify "Delta" Genes (Mutations or New Genes)
        delta_seq_map = {}
        final_energies = {}

        for gene, spec_seq in specific_seq_map.items():
            if gene == target_gene: continue

            # Check Cache
            if gene in cached_data:
                gen_energy, gen_seq = cached_data[gene]
                if spec_seq == gen_seq:
                    # MATCH: Reuse cached energy
                    final_energies[gene] = gen_energy
                    continue

            # MISMATCH: Sequence differs (or new gene) -> Add to Delta
            delta_seq_map[gene] = spec_seq

        # 4. Run RIsearch ONLY on Delta Genes
        if delta_seq_map:
            trigger = dna_to_rna_reverse_complement(aso_seq)
            delta_result = get_trigger_mfe_scores_by_risearch(
                trigger,
                delta_seq_map,
                minimum_score=cutoff,
                parsing_type='2'
            )
            delta_df = parse_risearch_output(delta_result)
            delta_agg = aggregate_off_targets(delta_df)

            # Add delta energies to final map
            if not delta_agg.empty:
                for _, r in delta_agg.iterrows():
                    final_energies[r['target']] = r['energy']

        # 5. Calculate Final Score
        score = calculate_score_helper(final_energies, specific_exp_map, method)
        index_scores[index] = score

    return index_scores


def calculate_score_helper(energy_dict, expression_dict, method):
    """
    Standardizes the scoring math to avoid code duplication.
    """
    score = 0.0
    RT = 0.616

    if not energy_dict:
        return 0.0

    # Pre-filter for valid interactions (Energy < 0) and Expression > 0
    # This speeds up the loops significantly
    valid_targets = []
    for gene, energy in energy_dict.items():
        exp = expression_dict.get(gene, 0)
        # Use 1e6 scaling for TPM-based methods (1, 2, 3, 4, 5)
        # Method 0 uses norm, which shouldn't be scaled by 1e6
        if method != 0:
            exp = exp / 1e6

        if energy < 0 and exp > 0:
            valid_targets.append((gene, energy, exp))

    if not valid_targets:
        return 0.0

    if method == 0:  # Normalized Arithmetic (Direct exp, no 1e6 scaling)
        # Re-loop because we scaled exp above, but method 0 usually takes raw norm
        score = 0.0
        for gene, energy in energy_dict.items():
            exp = expression_dict.get(gene, 0)  # Raw Norm
            if energy < 0 and exp > 0:
                score += energy * exp

    elif method == 1:  # TPM Arithmetic
        for _, energy, exp in valid_targets:
            score += energy * exp

    elif method == 2:  # Weighted Arithmetic (e^2)
        for _, energy, exp in valid_targets:
            score += energy * (exp ** 2)

    elif method == 3:  # Geometric
        log_sum = 0.0
        for _, energy, exp in valid_targets:
            log_sum += exp * np.log(-energy)
        score = log_sum

    elif method == 5:  # Mechanical Statistics
        for _, energy, exp in valid_targets:
            score += exp * np.exp(-energy / RT)

    elif method == 4:  # Ranked
        # Sort by expression descending
        valid_targets.sort(key=lambda x: x[2], reverse=True)

        for rank, (_, energy, exp) in enumerate(valid_targets, start=1):
            batch = (rank - 1) // 10
            weight = 10 / (2 ** batch)
            score += energy * exp * weight

    # Normalize Mechanical Statistics if required (usually done outside, but can be done here if context allows)
    # Note: Normalization requires knowing the max of the whole vector, so strictly speaking
    # method 5 normalization should happen after collecting all scores.
    # For now, this returns the raw Sum.

    return score