from asodesigner.fold import get_trigger_mfe_scores_by_risearch
from scripts.Roni.off_target_functions import dna_to_rna_reverse_complement, parse_risearch_output, aggregate_off_targets
import pandas as pd
import numpy as np
from io import StringIO
from asodesigner.consts import DATA_PATH_NEW
import os
import pickle

print('run')

# Set base directory relative to script
script_dir = os.path.dirname(__file__)
data_dir = os.path.abspath(os.path.join(script_dir, "..", "data_genertion"))

# Load the main ASO dataset
data_path = os.path.join(data_dir, "data_asoptimizer_updated.csv")

may_df = pd.read_csv(data_path)
#may_df = may_df.head(5)

print('read ASO data')

# Expression files path
expr_path = os.path.join(data_dir, "cell_line_expression")

# Load each transcriptome file -================== MATURE mRNA ===========================
'''
A431_df = pd.read_csv(os.path.join(expr_path, 'ACH-001328_transcriptome.csv'))
NCI_H460_df = pd.read_csv(os.path.join(expr_path, 'ACH-000463_transcriptome.csv'))
SH_SY5Y_df = pd.read_csv(os.path.join(expr_path, 'ACH-001188_transcriptome.csv'))
HeLa_df = pd.read_csv(os.path.join(expr_path, 'ACH-001086_transcriptome.csv'))
HepG2_df = pd.read_csv(os.path.join(expr_path, 'ACH-000739_transcriptome.csv'))
U_251MG_df = pd.read_csv(os.path.join(expr_path, 'ACH-000232_transcriptome.csv'))
'''
# Load each transcriptome file -================== PRE mRNA ===========================
# A431_df = pd.read_csv(DATA_PATH_NEW/'cell_line_expression/ACH-001328.mutated_transcriptome_premRNA.merged.csv')
# NCI_H460_df = pd.read_csv(DATA_PATH_NEW/'cell_line_expression/ACH-000463.mutated_transcriptome_premRNA.merged.csv')
# SH_SY5Y_df = pd.read_csv(DATA_PATH_NEW/'cell_line_expression/ACH-001188.mutated_transcriptome_premRNA.merged.csv')
# HeLa_df = pd.read_csv(DATA_PATH_NEW/'cell_line_expression/ACH-001086.mutated_transcriptome_premRNA.merged.csv')
# HepG2_df = pd.read_csv(DATA_PATH_NEW/'cell_line_expression/ACH-000739.mutated_transcriptome_premRNA.merged.csv')
# U_251MG_df = pd.read_csv(DATA_PATH_NEW/'cell_line_expression/ACH-000232.mutated_transcriptome_premRNA.merged.csv')



# ============================ Cut to top n expressed ==============================


# A431_df = A431_df.head(n)
# NCI_H460_df = NCI_H460_df.head(n)
# SH_SY5Y_df = SH_SY5Y_df.head(n)
# HeLa_df = HeLa_df.head(n)
# HepG2_df = HepG2_df.head(n)
# U_251MG_df = U_251MG_df.head(n)

 # ===================================================================================
cell_line_list = ['ACH-001328', 'ACH-000463', 'ACH-001188', 'ACH-001086', 'ACH-000739', 'ACH-000232']

# cell_line2df_ = {
#     'A431':A431_df,
#     'A-431':A431_df,
#     'NCI-H460':NCI_H460_df,
#     'SH_SY5Y':SH_SY5Y_df,
#     'SH-SY5Y':SH_SY5Y_df,
#     'HeLa':HeLa_df,
#     'Hela':HeLa_df,
#     'HepG2':HepG2_df,
#     'U-251MG':U_251MG_df,
#     'U251':U_251MG_df
# }

cellline_to_filename = {
    'A431': 'ACH-001328',
    'A-431': 'ACH-001328',
    'NCI-H460': 'ACH-000463',
    'SH_SY5Y': 'ACH-001188',
    'SH-SY5Y': 'ACH-001188',
    'HeLa': 'ACH-001086',
    'Hela': 'ACH-001086',
    'HepG2': 'ACH-000739',
    'U-251MG': 'ACH-000232',
    'U251': 'ACH-000232'
}


def get_off_target_info(cellline_to_filename, ASO_df, cutoff, topN):
    index_info_vec = {}

    # Loop over each cell line once
    for cell_line, group_df in ASO_df.groupby("Cell_line"):

        # skip if cell line not in map
        if cell_line not in cellline_to_filename:
            print(f"Skipping {cell_line}, not in mapping")
            continue

        # load transcriptome data once for this cell line
        cell_line_id = cellline_to_filename[cell_line]
        transcriptome_path = DATA_PATH_NEW / f"cell_line_expression/{cell_line_id}.mutated_transcriptome_premRNA.merged.csv"
        if not transcriptome_path.exists():
            print(f"File missing for {cell_line}: {transcriptome_path}")
            continue

        curr_df = pd.read_csv(transcriptome_path).head(topN)
        print(f"Loaded {cell_line} ({cell_line_id}), {len(curr_df)} transcripts")

        # preprocess transcriptome into dictionaries
        name_to_seq = {}
        name_to_exp_TPM = {}
        name_to_exp_norm = {}

        for _, gene_row in curr_df.iterrows():
            curr_gene = gene_row['Gene'].split()[0]

            mut_seq = gene_row.get('Mutated Transcript Sequence')
            og_seq = gene_row.get('Original Transcript Sequence')
            if pd.isna(mut_seq) and pd.isna(og_seq):
                continue

            mRNA_seq = mut_seq if not pd.isna(mut_seq) else og_seq
            name_to_seq[curr_gene] = mRNA_seq
            name_to_exp_TPM[curr_gene] = gene_row['expression_TPM']
            name_to_exp_norm[curr_gene] = gene_row['expression_norm']

        # now loop over all ASOs in this cell line
        for idx, row in group_df.iterrows():
            index = row['index']
            ASO_seq = row['Sequence']
            trigger = dna_to_rna_reverse_complement(ASO_seq)
            target_gene = row['Canonical Gene Name']

            # drop the target gene from search space
            seqs = {g: s for g, s in name_to_seq.items() if g != target_gene}
            exp_TPM = {g: e for g, e in name_to_exp_TPM.items() if g != target_gene}
            exp_norm = {g: e for g, e in name_to_exp_norm.items() if g != target_gene}

            # run mfe search
            result_dict = get_trigger_mfe_scores_by_risearch(
                trigger,
                seqs,
                minimum_score=cutoff,
                parsing_type='2'
            )

            result_df = parse_risearch_output(result_dict)

            if result_df.empty:
                index_info_vec[index] = pd.DataFrame()
                continue

            # keep best score per target (might change later to MECHANICAL STATISTICS-based binding)
            top_score_df = result_df.loc[result_df.groupby('target')['score'].idxmax()]

            # add weighted expression features
            g = top_score_df['energy']
            e = top_score_df['target'].map(exp_TPM)/1e6

            top_score_df['energy_w_by_exp_ARTH'] = g * e
            top_score_df['energy_w_by_log_exp_ARTH'] = g * top_score_df['target'].map(exp_norm)
            top_score_df['energy_w_by_exp_GEO'] = -((-g)**e)
            top_score_df['energy_w_by_log_exp_GEO'] = - e * np.log(-g)

            top_score_df = top_score_df.sort_values(by='energy_w_by_log_exp_GEO',
                                                    ascending=True).reset_index(drop=True)

            index_info_vec[index] = top_score_df

    return index_info_vec

n = 100
cutoff = 800

off_target_info_vec = get_off_target_info(cellline_to_filename, may_df, cutoff, n)
with open(f'off_target_info.premRNA.top{n}.cutoff{cutoff}.pkl', 'wb') as f:
    pickle.dump(off_target_info_vec, f)
print("yay saved")
# off_target_feature = pd.DataFrame({
#     "off_target_score_TPM": off_target_vec[0],
# }).reset_index().rename(columns={"index": "index"})

#off_target_feature.to_csv(os.path.join(data_dir, "off_target___.csv"), index=False)
#print("off_target_feature.csv saved")

