import numpy as np
from numba import njit

# NOTE: This code matches exactly the codonbias implementation, but is much quicker asof 0.3.5 version


@njit
def char_to_int(c: int) -> int:
    if c == 65 or c == 97:
        return 0  # A
    if c == 67 or c == 99:
        return 1  # C
    if c == 71 or c == 103:
        return 2  # G
    if c == 84 or c == 116 or c == 85 or c == 117:
        return 3  # T/U
    return -1


# Mappings equivalent to the dictionaries above (NCBI ID 1)
CODON_TO_AA = np.array(
    [
        8,
        11,
        8,
        11,
        16,
        16,
        16,
        16,
        14,
        15,
        14,
        15,
        7,
        7,
        10,
        7,
        13,
        6,
        13,
        6,
        12,
        12,
        12,
        12,
        14,
        14,
        14,
        14,
        9,
        9,
        9,
        9,
        3,
        2,
        3,
        2,
        0,
        0,
        0,
        0,
        5,
        5,
        5,
        5,
        17,
        17,
        17,
        17,
        -1,
        19,
        -1,
        19,
        15,
        15,
        15,
        15,
        -1,
        1,
        18,
        1,
        9,
        4,
        9,
        4,
    ],
    dtype=np.int32,
)

AA_DEG = np.array(
    [4, 2, 2, 2, 2, 4, 2, 3, 2, 6, 1, 2, 4, 2, 6, 6, 4, 4, 1, 2], dtype=np.int32
)
DEG_VALUES = np.array([1, 2, 3, 4, 6], dtype=np.int32)
DEG_COUNTS = np.array([2, 9, 1, 5, 3], dtype=np.float64)


@njit
def _compute_enc_core(seq_bytes) -> float:
    n_len = len(seq_bytes)
    if n_len == 0:
        return np.nan

    # 1. Background Nucleotide Composition
    # CRITICAL FIX: start with 1.0 pseudocount to match BaseCounter.get_table(normed=True)
    bnc_counts = np.ones(4, dtype=np.float64)
    cod_counts = np.zeros(64, dtype=np.float64)

    for i in range(n_len):
        b = char_to_int(seq_bytes[i])
        if b >= 0:
            bnc_counts[b] += 1.0

    bnc_sum = np.sum(bnc_counts)
    BNC = bnc_counts / bnc_sum

    for i in range(0, n_len - 2, 3):
        b1 = char_to_int(seq_bytes[i])
        b2 = char_to_int(seq_bytes[i + 1])
        b3 = char_to_int(seq_bytes[i + 2])
        if b1 >= 0 and b2 >= 0 and b3 >= 0:
            idx = b1 * 16 + b2 * 4 + b3
            if CODON_TO_AA[idx] != -1:
                cod_counts[idx] += 1.0

    # 2. Background Codon Composition (BCC)
    BCC = np.zeros(64, dtype=np.float64)
    for idx in range(64):
        aa = CODON_TO_AA[idx]
        if aa != -1:
            b1 = idx // 16
            b2 = (idx // 4) % 4
            b3 = idx % 4
            BCC[idx] = BNC[b1] * BNC[b2] * BNC[b3]

    for a in range(20):
        sum_bcc = 0.0
        for idx in range(64):
            if CODON_TO_AA[idx] == a:
                sum_bcc += BCC[idx]
        if sum_bcc > 0.0:
            for idx in range(64):
                if CODON_TO_AA[idx] == a:
                    BCC[idx] /= sum_bcc

    # 3. Add pseudocount and calculate Chi-Squared (F)
    cod_counts += 1.0
    F_aa = np.zeros(20, dtype=np.float64)
    N_aa = np.zeros(20, dtype=np.float64)

    for a in range(20):
        n_sum = 0.0
        for idx in range(64):
            if CODON_TO_AA[idx] == a:
                n_sum += cod_counts[idx]
        N_aa[a] = n_sum

        chi2_a = 0.0
        deg = AA_DEG[a]
        for idx in range(64):
            if CODON_TO_AA[idx] == a:
                P_idx = cod_counts[idx] / n_sum
                diff = P_idx - BCC[idx]
                chi2_a += (diff * diff) / BCC[idx]

        F_aa[a] = (chi2_a + 1.0) / deg

    # 4. Weighted mean of F across degeneracy groups
    F_deg_sum = np.zeros(5, dtype=np.float64)
    N_deg_sum = np.zeros(5, dtype=np.float64)

    for a in range(20):
        deg = AA_DEG[a]
        deg_idx = (
            0
            if deg == 1
            else 1
            if deg == 2
            else 2
            if deg == 3
            else 3
            if deg == 4
            else 4
        )
        F_deg_sum[deg_idx] += F_aa[a] * N_aa[a]
        N_deg_sum[deg_idx] += N_aa[a]

    F_deg = np.zeros(5, dtype=np.float64)
    for i in range(5):
        if N_deg_sum[i] > 0.0:
            F_deg[i] = F_deg_sum[i] / N_deg_sum[i]
        else:
            F_deg[i] = np.nan

    if np.isnan(F_deg[2]):  # Fix missing deg=3
        F_deg[2] = 0.5 * (F_deg[1] + F_deg[3])

    # 5. Final Calculation
    ENC = 0.0
    for i in range(5):
        ENC += DEG_COUNTS[i] / F_deg[i]

    return 61.0 if ENC > 61.0 else ENC


# TODO: merge this logic to codo-bias and then return to the simple function call


def compute_ENC(seq: str) -> float:
    if not isinstance(seq, str) or not seq.strip():
        return np.nan

    seq_bytes = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    enc_score = _compute_enc_core(seq_bytes)

    if np.isnan(enc_score):
        return np.nan

    return (enc_score - 20.0) / (61.0 - 20.0)
