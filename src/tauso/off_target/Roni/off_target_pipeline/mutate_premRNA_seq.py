import pandas as pd
import pickle
from mutate_cell_line_transcriptome import mutate, mutation_dict, celline_list, find_shift

print('test')
"""
I've made this code because when I switched to working on pre-mRNA 
rather than mature mRNA the final_func didn't run due to RAM limitations.
"""


path = '/home/oni/ASOdesign/scripts/data_genertion/cell_line_expression/'
addition_exp = '_transcriptome_premRNA.merged.csv'

with open(path + '_gtf_annotations.pkl', 'rb') as f:
    annotations = pickle.load(f)
print('GTF annotations loaded.')

for cell_line in celline_list:
    print(f'processing {cell_line}')
    # Extract necessary data
    mut_data = pd.read_csv(path + cell_line + '_mutations.csv')
    exp_data = pd.read_csv(path + cell_line + addition_exp)

    exp_data_indexed = exp_data.set_index('Transcript_ID')
    print('Data loaded successfully.')

    for idx, row in mut_data.iterrows():
        mut_dict = mutation_dict(row)

        if mut_dict is not None:
            tid = mut_dict['id']
            mut_idx = mut_dict['start']


            if tid in exp_data_indexed.index:
                try:
                    shift = find_shift(annotations[tid], mut_idx)
                    print(f"{tid} shift: {shift}")
                    if pd.isna(exp_data_indexed.at[tid, 'Mutated Transcript Sequence']):
                        seq = exp_data_indexed.at[tid, 'Original Transcript Sequence']

                    else:
                        seq = exp_data_indexed.at[tid, 'Mutated Transcript Sequence']

                    mutated_seq = mutate(mut_dict, shift, seq)
                    exp_data_indexed.at[tid, 'Mutated Transcript Sequence'] = mutated_seq
                    if seq != mutated_seq:
                        print(f'really mutated {tid}')
                    else:
                        print(f'didnt really do anything... check {tid}')
                except Exception as e:
                    print(f"Skipping {tid} due to error: {e}")

    exp_data_indexed.to_csv(path + cell_line + '.mutated' + addition_exp)
    print(f'Saved {cell_line}')
