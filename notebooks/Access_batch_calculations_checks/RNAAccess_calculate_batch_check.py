import time
import numpy as np
import sys
import os
from pathlib import Path
current_file_path = Path(__file__).resolve()
project_root = current_file_path.parents[2]
sys.path.append(str(project_root))
from tauso.features.rna_access.rna_access import RNAAccess
from tauso._raccess.core import find_raccess


YEAST_M_CHERRY = (
    "ATGTCTAAGGGGGAAGAAGACAATATGGCGATTATTAAAGAGTTTATGAGATTTAAAGTACATATGGAAGGAAGTGTTAAT"
    "GGTCACGAGTTTGAGATCGAAGGTGAAGGTGAAGGTCGTCCATATGAGGGTACGCAAACAGCAAAACTAAAGGTGACTAAAG"
    "GGGGACCATTACCTTTCGCTTGGGATATACTGTCACCACAATTCATGTACGGATCGAAAGCTTACGTAAAGCACCCGGCCGA"
).replace("T", "U")

def build_two_seqs(gene):
    # building two sequences out of a gene for testing
    seq1 = gene[0:120]
    seq2 = gene[30:150]
    return [("rna1", seq1), ("rna2", seq2)]

def build_n_seqs(n,gene):
    length = len(gene)
    if length < n+120:
        return "too short gene - cant provide n seqs"
    seqs = []
    for i in range(n):
        seq = gene[i:120+i]
        seqs.append((f"rna{i}", seq))
    return seqs

def time_batch_vs_single(n):
    segment_sizes = [8]  # check if that's a relevant value
    max_span = 120

    exe_path = find_raccess()
    ra = RNAAccess(exe_path=exe_path, segment_sizes=segment_sizes, max_span=max_span)

    seq_id_list = build_n_seqs(n,YEAST_M_CHERRY)

    t0 = time.perf_counter()
    res_batch = ra.calculate(seq_id_list)
    t_batch = time.perf_counter() - t0
    print(f"time of {n} batch size is {t_batch:.3f} seconds")

    res_single = {}
    t0 = time.perf_counter()
    for pair in seq_id_list:
        tmp = ra.calculate([pair])
        res_single[pair[0]] = tmp[pair[0]]
    t_single = time.perf_counter() - t0
    print(f"{n} individual calls took {t_single:.3f} seconds")

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

def main():
    segment_sizes = [8] # check if that's a relevant value
    max_span = 120

    exe_path = find_raccess()
    ra = RNAAccess(exe_path = exe_path, segment_sizes = segment_sizes, max_span = max_span)

    seq_id_list = build_two_seqs(YEAST_M_CHERRY)

    print("running raaccess on 2 seqs")
    res = ra.calculate(seq_id_list)
    print(res)
    print("\n=== Timing batch vs single ===")
    time_batch_vs_single(8)

if __name__ == "__main__":
    main()