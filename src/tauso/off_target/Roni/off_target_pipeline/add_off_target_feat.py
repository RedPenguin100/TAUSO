from src.tauso.hybridization.fast_hybridization import get_trigger_mfe_scores_by_risearch
from src.tauso.off_target.Roni.off_target_pipeline.off_target_functions import dna_to_rna_reverse_complement, parse_risearch_output, aggregate_off_targets

import pandas as pd
import numpy as np



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


def add_off_target_feat(cell_line_dict, cell_line2df, ASO_df, cutoff, method=0):
    index_score_vec_TPM = {}
    index_score_vec_norm = {}
    index_score_vec_MS = {}
    index_score_vec_geo = {}
    index_score_vec = {}

    for idx, row in ASO_df.iterrows(): # For ASO

        index = row['index']
        exp_cell_line = row['Cell_line']
        cell_line = cell_line_dict[exp_cell_line]
        ASO_seq = row['Sequence']
        trigger = dna_to_rna_reverse_complement(ASO_seq)
        target_gene = row['Canonical Gene Name']

        if cell_line not in cell_line2df:
            continue

        curr_df = cell_line2df[cell_line]
        print(curr_df.head())

        name_to_seq = {}
        name_to_exp_TPM = {}
        name_to_exp_norm = {}

        for _, gene_row in curr_df.iterrows(): # For transcript
            curr_gene = gene_row['Gene'].split()[0]

            # Skip the target gene
            if curr_gene == target_gene:
                continue

            # Select mutated or original sequence
            mut_seq = gene_row.get('Mutated Transcript Sequence')
            og_seq = gene_row.get('Original Transcript Sequence')

            if pd.isna(mut_seq) and pd.isna(og_seq):
                continue

            mRNA_seq = mut_seq if not pd.isna(mut_seq) else og_seq
            exp_TPM = gene_row['expression_TPM']
            exp_norm = gene_row[f'{cell_line}_expression_norm']

            name_to_seq[curr_gene] = mRNA_seq
            name_to_exp_TPM[curr_gene] = exp_TPM
            name_to_exp_norm[curr_gene] = exp_norm

        if not name_to_seq:
            print(f"[WARNING] name_to_seq is empty for index={index} cell_line={cell_line}")

        # Calculate mfe scores per ASO to all topN expressed genes
        result_dict = get_trigger_mfe_scores_by_risearch(
            trigger,
            name_to_seq,
            minimum_score=cutoff,
            parsing_type='2'
        )
        result_df = parse_risearch_output(result_dict) # make it interperable
        #print(f'001:\n {result_df}')
        result_df_agg = aggregate_off_targets(result_df) # take the maximal bind per transcript and make a proper df
        #print(f'002:\n {result_df_agg}')

        # ===================== Weight by expression (Arithmetic mean) ==========================
        if method == 0:
            result_dict_weighted_TPM = {
                row['target']: row['energy'] * name_to_exp_TPM.get(row['target'], 0)
                for _, row in result_df_agg.iterrows()
            }
            result_dict_weighted_norm = {
                row['target']: row['energy'] * name_to_exp_norm.get(row['target'], 0)
                for _, row in result_df_agg.iterrows()
            }
            # print(result_dict_weighted_TPM)
            index_score_vec_TPM[index] = sum(result_dict_weighted_TPM.values())

            # print(result_dict_weighted_norm)
            index_score_vec_norm[index] = sum(result_dict_weighted_norm.values())

        # ================== Weight by expression (Geometric mean) =============================
        if method == 1:
            log_sum = 0.0
            count = 0
            for _, r in result_df_agg.iterrows():
                g = r['energy']
                e = float(name_to_exp_TPM.get(r['target'], 0)) / 1e6
                if (g is not None) and (e > 0) and (g < 0):
                    log_sum += e * np.log(-g)
                    count += 1

            product_log = log_sum
            product = float(np.exp(log_sum)) if count > 0 else 0.0
            index_score_vec_geo[index] = product_log

    # =================== Weight by weighted expression arithmetic mean ===================================
        if method == 2:
            _sum = 0
            power = 2
            for _, r in result_df_agg.iterrows():
                g = r['energy']
                e = float(name_to_exp_TPM.get(r['target'], 0)) / 1e6
                if (g is not None) and (e > 0) and (g < 0):
                    _sum += g * e**power

            index_score_vec[index] = _sum

    # =================== Linear decreasing Arithmetic weighting (jumps of 10) ========================================
        if method == 3:
            # sort transcripts by expression (highest first)
            ranked = result_df_agg.assign(
                exp=[name_to_exp_TPM.get(t, 0) for t in result_df_agg['target']]
            ).sort_values("exp", ascending=False)

            _sum = 0
            for rank, (_, r) in enumerate(ranked.iterrows(), start=1):
                g = r['energy']
                e = float(r['exp']) / 1e6

                if (g is not None) and (e > 0) and (g < 0):
                    batch = (rank - 1) // 10  # 0 for top 1–10, 1 for 11–20, etc.
                    weight = 10 / (2 ** batch)
                    _sum += g * e * weight

            index_score_vec[index] = _sum


    # =================== Weight by mechanical statistics ===================================
        if method == 4:
            RT = 0.616
            bind_score = 0
            for _, r in result_df_agg.iterrows():
                g = r['energy']
                e = float(name_to_exp_TPM.get(r['target'], 0)) / 1e6
                bind_score += e * np.exp(-g / RT)

            index_score_vec_MS[index] = bind_score

    if method == 4:
        max_val = max(index_score_vec_MS.values())
        if max_val > 0:
            for k, v in index_score_vec_MS.items():
                index_score_vec_MS[k] = v / max_val
        return index_score_vec_MS

    #return [index_score_vec_TPM, index_score_vec_norm]
    return index_score_vec

'''
off_target_vec = get_off_target_feature(cell_line2df_, may_df)
#print(off_target_vec)

off_target_feature = pd.DataFrame({
    "off_target_score_TPM": off_target_vec[0],
    "off_target_score_log": off_target_vec[1],
}).reset_index().rename(columns={"index": "index"})

off_target_feature.to_csv(os.path.join(data_dir, "off_target.top500.cutoff1200.premRNA.csv"), index=False)
print("off_target_feature.csv saved")

'''