import pandas as pd
from mutate_cell_line_transcriptome import celline_list

path = '/home/oni/ASOdesign/scripts/data_genertion/cell_line_expression/'
sequences = '.mutated_transcriptome_premRNA.merged.csv'
muts = '_mutations.csv'

def compare_sequences(seq1, seq2):
    """
    Compare two strings and report where they differ.
    Indices are 0-based.
    """
    len1, len2 = len(seq1), len(seq2)
    min_len = min(len1, len2)

    # Find first difference from left
    left_diff = None
    for i in range(min_len):
        if seq1[i] != seq2[i]:
            left_diff = i
            break

    # Find first difference from right
    right_diff = None
    for i in range(1, min_len + 1):
        if seq1[-i] != seq2[-i]:
            right_diff = max(len1, len2) - i
            break

    # If sequences are identical
    if left_diff is None and len1 == len2:
        print("Sequences are identical.")
        return

    # Slice out the differing parts
    diff_seq1 = seq1[left_diff:right_diff+1] if left_diff is not None else ""
    diff_seq2 = seq2[left_diff:right_diff+1] if left_diff is not None else ""
    area_seq1 = seq1[left_diff-30:right_diff+30] if left_diff is not None else ""
    area_seq2 = seq2[left_diff-30:right_diff+30] if left_diff is not None else ""

    print(f"First difference at index: {left_diff}")
    print(f"Last difference at index: {right_diff}")
    print(f"Seq1 diff: {diff_seq1}")
    print(f"Seq2 diff: {diff_seq2}")
    print(f"Seq1 area: {area_seq1}")
    print(f"Seq2 area: {area_seq2}\n")

def do_da_compare():
    for cellline in celline_list:
        seq_data, mut_data = pd.read_csv(path + cellline + sequences), pd.read_csv(path + cellline + muts)
        print(f"read {cellline} data")
        for idx, row in mut_data.iterrows():
            tid, mut = row['DNAChange'].split(":")
            for i, r in seq_data.iterrows():
                if r['Transcript_ID'] == tid:
                    print(f"for {tid} and change {mut}:")
                    compare_sequences(r['Original Transcript Sequence'], r['Mutated Transcript Sequence'])

def turn_T_U(seq):
    return seq.replace('T', 'U')

for cellline in celline_list:
    seq_data = pd.read_csv(path + cellline + sequences)
    print(f"Read {cellline} data")

    seq_data['Original Transcript Sequence'] = seq_data['Original Transcript Sequence'].apply(turn_T_U)
    seq_data['Mutated Transcript Sequence'] = seq_data['Mutated Transcript Sequence'].apply(
        lambda x: turn_T_U(x) if not pd.isna(x) else x
    )

    seq_data.to_csv(path + cellline + sequences, index=False)