import time
import numpy as np
import pandas as pd
import sys
import os
from pathlib import Path
current_file_path = Path(__file__).resolve()
project_root = current_file_path.parents[2]
sys.path.append(str(project_root))
from notebooks.consts import *
from notebooks.notebook_utils import log_correction, read_cached_gene_to_data
from tauso.util import get_antisense
from tauso.features.rna_access.access_calculator import AccessCalculator
from RNAAccess_calculate_batch_check import build_n_seqs
from tauso.features.rna_access.rna_access import RNAAccess
from tauso.features.rna_access.access_calculator import get_sense_with_flanks



YEAST_M_CHERRY = (
    "ATGTCTAAGGGGGAAGAAGACAATATGGCGATTATTAAAGAGTTTATGAGATTTAAAGTACATATGGAAGGAAGTGTTAAT"
    "GGTCACGAGTTTGAGATCGAAGGTGAAGGTGAAGGTCGTCCATATGAGGGTACGCAAACAGCAAAACTAAAGGTGACTAAAG"
    "GGGGACCATTACCTTTCGCTTGGGATATACTGTCACCACAATTCATGTACGGATCGAAAGCTTACGTAAAGCACCCGGCCGA"
).replace("T", "U")

def check_batch_vs_single_calc(n):
    seed_sizes = [8]
    max_span = 120
    access_size = 13

    seq_id_list = build_n_seqs(n, YEAST_M_CHERRY)

    t0 = time.perf_counter()
    res_batch = AccessCalculator.calc_access_energies_batch(
            seq_id_list, access_size, seed_sizes, max_span
    )
    t1 = time.perf_counter()


    res_single = {}
    t2 = time.perf_counter()
    for rna_id, seq in seq_id_list:
        df = AccessCalculator.calc_access_energies(
            seq, access_size, seed_sizes, max_span, uuid_str = None)
        res_single[rna_id] = df
    t3 = time.perf_counter()

    for rna_id in res_batch.keys():
        df_b = res_batch[rna_id]
        df_s = res_single[rna_id]

        # shape, index, columns identical?
        assert df_b.shape == df_s.shape
        assert (df_b.index == df_s.index).all()
        assert (df_b.columns == df_s.columns).all()

        # numeric values identical within floating-point tolerance
        assert np.allclose(df_b.values, df_s.values, rtol=1e-8, atol=1e-10)

        print("All batch results match single results exactly!")

    print(f'batch calculation took {t1 - t0} seconds')
    print(f'single calculation took {t3 - t2} seconds')

csv_path = Path("..") / "data" / "data_asoptimizer_updated.csv"
all_data = pd.read_csv(str(csv_path), low_memory=False)
# Remove rows with missing values in the INHIBITION column
all_data_no_nan = all_data.dropna(subset=[INHIBITION]).copy()
# Create a new column with transformed inhibition values on a negative log scale
log_correction(all_data_no_nan) # to avoid log 0
# Filter the data to include only rows where the cell line organism is human
all_data_no_nan_human = all_data_no_nan[all_data_no_nan[CELL_LINE_ORGANISM] == 'human']
genes = all_data_no_nan[CANONICAL_GENE].copy()
genes_u = list(set(genes))
# Remove non-human or negative controls from the gene list
genes_u.remove('HBV')
genes_u.remove('negative_control')
gene_to_data = read_cached_gene_to_data(genes_u)
# Filter data to keep only rows with valid gene information
all_data_human_gene = all_data_no_nan_human[all_data_no_nan_human[CANONICAL_GENE].isin(genes_u)].copy()

# Define names for new columns
SENSE_SEQUENCE = 'sense_sequence'
PRE_MRNA_SEQUENCE = 'pre_mrna_sequence'
SENSE_START = 'sense_start'
SENSE_LENGTH = 'sense_length'

# Initialize new columns
all_data_human_gene[SENSE_SEQUENCE] = ""
all_data_human_gene[PRE_MRNA_SEQUENCE] = ""
all_data_human_gene[SENSE_START] = np.zeros_like(all_data_human_gene[CANONICAL_GENE], dtype=int)
all_data_human_gene[SENSE_LENGTH] = np.zeros_like(all_data_human_gene[CANONICAL_GENE], dtype=int)

# Iterate over each row and compute the antisense complement and the gene's pre-mRNA
for index, row in all_data_human_gene.iterrows():
    gene_name = row[CANONICAL_GENE]

    if gene_name not in gene_to_data:
        continue  # Skip genes not found in genome annotation

    locus_info = gene_to_data[gene_name]
    pre_mrna = locus_info.full_mrna
    antisense = row[SEQUENCE]
    sense = get_antisense(antisense)
    idx = pre_mrna.find(sense)

    # Store computed sequences in new columns
    all_data_human_gene.loc[index, SENSE_START] = idx
    all_data_human_gene.loc[index, SENSE_LENGTH] = len(antisense)
    all_data_human_gene.at[index, SENSE_SEQUENCE] = sense
    all_data_human_gene.at[index, PRE_MRNA_SEQUENCE] = pre_mrna

FLANK_SIZE = 120  # Change this as needed
FLANKED_SENSE_COL = f'sense_with_flank_{FLANK_SIZE}nt'

# Create new column with flanked sequences
all_data_human_gene[FLANKED_SENSE_COL] = all_data_human_gene.apply(
    lambda row: get_sense_with_flanks(
        row['pre_mrna_sequence'],
        row['sense_start'],
        row['sense_length'],
        flank_size=FLANK_SIZE
    ) if row['sense_start'] != -1 else "",  # Handle cases where sense was not found
    axis=1
)
valid_data = all_data_human_gene[
    all_data_human_gene[FLANKED_SENSE_COL].astype(str).str.len() > 0
].copy()

def build_flank_seq_list(df, n, flank_col=FLANKED_SENSE_COL):
    df_sample = df.sample(n=n, random_state=42)
    seq_id_list = []
    assert "index" in df.columns
    for _, row in df_sample.iterrows():
        seq = row[flank_col]
        seq_id = str(row["index"])
        seq_id_list.append((seq_id,seq.upper().replace('T','U')))
    return df_sample, seq_id_list

def check_batch_vs_single_calc_real_data(df,n):
    seed_sizes = [8]
    max_span = 120
    access_size = 13

    _, seq_id_list = build_flank_seq_list(df, n, flank_col=FLANKED_SENSE_COL)

    t0 = time.perf_counter()
    res_batch = AccessCalculator.calc_access_energies_batch(
            seq_id_list, access_size, seed_sizes, max_span
    )
    t1 = time.perf_counter()


    res_single = {}
    t2 = time.perf_counter()
    for rna_id, seq in seq_id_list:
        df = AccessCalculator.calc_access_energies(
            seq, access_size, seed_sizes, max_span, uuid_str = None)
        res_single[rna_id] = df
    t3 = time.perf_counter()

    for rna_id in res_batch.keys():
        df_b = res_batch[rna_id]
        df_s = res_single[rna_id]

        # shape, index, columns identical?
        assert df_b.shape == df_s.shape
        assert (df_b.index == df_s.index).all()
        assert (df_b.columns == df_s.columns).all()

        # numeric values identical within floating-point tolerance
        assert np.allclose(df_b.values, df_s.values, rtol=1e-8, atol=1e-10)

    print("All batch results in sampled data match single results exactly!")

    print(f'batch calculation took for sampled data {t1 - t0} seconds')
    print(f'single calculation took sampled data {t3 - t2} seconds')


if __name__ == "__main__":
    check_batch_vs_single_calc(10)
    check_batch_vs_single_calc_real_data(valid_data, 50)
