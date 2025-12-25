import pandas as pd
import os
from glob import glob
# from get_premRNA_sequences_II import cell_line2dataII



path_0 = '/home/oni/ASOdesign/scripts/data_genertion/cell_line_expression/'
path = '/home/oni/ASOdesign/scripts/Roni/'
cell_line_list = ['ACH-000232', 'ACH-000463', 'ACH-000739', 'ACH-001086', 'ACH-001188', 'ACH-001328', 'ACH-000681']
addition = '_transcriptome_premRNA'



for n in range(6, 7):
    print(n)
    df1 = pd.read_csv(path_0 + cell_line_list[n] + addition + '.csv')
    print(df1.shape)
    df2 = pd.read_csv(path + cell_line_list[n] + addition + '.chr_5_8.csv')
    df3 = pd.read_csv(path + cell_line_list[n] + addition + '.chr_9_12.csv')
    df4 = pd.read_csv(path + cell_line_list[n] + addition + '.chr_13_16.csv')
    df5 = pd.read_csv(path + cell_line_list[n] + addition + '.chr_17_20.csv')
    df6 = pd.read_csv(path + cell_line_list[n] + addition + '.chr_21_22_XY.csv')
    df7 = pd.read_csv(path + cell_line_list[n] + addition + '.chr_M_rest.csv')

    df_list = [df2, df3, df4, df5, df6, df7]  # exclude df1 itself

    df1.set_index("Transcript_ID", inplace=True)
    for df in df_list:
        print(df.shape)
        df.set_index("Transcript_ID", inplace=True)
        for col in df.columns:
            if col in df1.columns:
                df1[col] = df1[col].combine_first(df[col])
            else:
                df1[col] = df[col]
        print('merged')

    df1.reset_index(inplace=True)
    df1_sorted = df1.sort_values(by=['expression_TPM'], ascending=False)
    print(df1_sorted.shape)
    df1_sorted.to_csv(path + cell_line_list[n] + addition + '.merged.csv', index=False)

def merge_cell_line(df_main, df_lst):
    df_main.set_index("Transcript_ID", inplace=True)
    for df in df_lst:
        df.set_index("Transcript_ID", inplace=True)
        for col in df.columns:
            if col in df_main.columns:
                df_main[col] = df_main[col].combine_first(df[col])
            else:
                df_main[col] = df[col]
        print('merged')
