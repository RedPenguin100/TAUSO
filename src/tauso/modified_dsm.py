import numpy as np
from external.risearch.RIsearch1.numpy_to_csv import dsm_variable_to_csv, numpy_to_csv

position_dict = {
    "AA": 0,
    "AC": 1,
    "AG": 2,
    "AU": 3,
    "AN": 4,
    "A-": 5,
    "CA": 6,
    "CC": 7,
    "CG": 8,
    "CU": 9,
    "CN": 10,
    "C-": 11,
    "GA": 12,
    "GC": 13,
    "GG": 14,
    "GU": 15,
    "GN": 16,
    "G-": 17,
    "UA": 18,
    "UC": 19,
    "UG": 20,
    "UU": 21,
    "UN": 22,
    "U-": 23,
    "NA": 24,
    "NC": 25,
    "NG": 26,
    "NU": 27,
    "NN": 28,
    "N-": 29,
    "-A": 30,
    "-C": 31,
    "-G": 32,
    "-U": 33,
    "-N": 34,
    "--": 35
}


def modify_single_weight(weight_matrix, modification):
    left_nc, right_nc, weight = modification
    reshaped = weight_matrix.reshape((36, 36))
    row = position_dict[left_nc]
    col = position_dict[right_nc]
    reshaped[row][col] = weight


# got error when did option 2 np_array.transpose(3, 2, 1, 0)
def rna_dna_to_dna_rna(np_array):
    # print(np_array[1][2][3][4])
    a = np_array.transpose(2, 3, 0, 1)
    # print(a[1][2][3][4])
    return a


def rna_dna_to_ps_dna_rna(np_array, eps_array):
    return rna_dna_to_dna_rna(np_array) - eps_array


def make_dsm_dna_rna(dsm_base='dsm_su95_rev_woGU_pos'):
    dsm_variable_to_csv(dsm_name=dsm_base, func=rna_dna_to_dna_rna, kwargs=None,
                        func_extend=rna_dna_to_dna_rna, kwargs_extend=None)


def make_dsm_ps_dna_rna(dsm="dsm_su95_rev_woGU_pos", eps_array=None, base=20):
    if eps_array is None:
        eps_array = np.full((6, 6, 6, 6), fill_value=base)
    dsm_variable_to_csv(dsm_name=dsm, func=rna_dna_to_ps_dna_rna, kwargs={'eps_array': -eps_array},
                        func_extend=rna_dna_to_dna_rna, kwargs_extend=None)

def mini_transpose(eps_array):
    reshaped = eps_array.reshape((36, 36))
    copy = reshaped.copy()
    nucleotides = ['A', 'C', 'G', 'U']
    for N1 in nucleotides:
        for N2 in nucleotides:
            base = position_dict[N1 + N2]
            target = position_dict[N2 + N1]
            reshaped[base][target] = copy[target][base] # copy necessary because we override


def experimental_weights_50_deg(eps_array, neg=False):
    if not neg:
        modify_single_weight(eps_array, ('AA', 'UU', 31))
        modify_single_weight(eps_array, ('AU', 'UA', 14))
        modify_single_weight(eps_array, ('AC', 'UG', 6))
        modify_single_weight(eps_array, ('AG', 'UC', 23))
        modify_single_weight(eps_array, ('UA', 'AU', 28))
        modify_single_weight(eps_array, ('UU', 'AA', 34))
        modify_single_weight(eps_array, ('UC', 'AG', 24))
        modify_single_weight(eps_array, ('UG', 'AC', 18))
        modify_single_weight(eps_array, ('CA', 'GU', 13))
        modify_single_weight(eps_array, ('CU', 'GA', 28))
        modify_single_weight(eps_array, ('CC', 'GG', 18))
        modify_single_weight(eps_array, ('CG', 'GC', 41))
        modify_single_weight(eps_array, ('GA', 'CU', 3))
        modify_single_weight(eps_array, ('GU', 'CA', 33))
        modify_single_weight(eps_array, ('GC', 'CG', 61))
        modify_single_weight(eps_array, ('GG', 'CC', 0))

        modify_single_weight(eps_array, ('-A', '-U', -106))
        modify_single_weight(eps_array, ('A-', 'U-', -102))
        modify_single_weight(eps_array, ('-U', '-A', -107))
        modify_single_weight(eps_array, ('U-', 'A-', -105))
        modify_single_weight(eps_array, ('-C', '-G', -121))
        modify_single_weight(eps_array, ('C-', 'G-', -98))
        modify_single_weight(eps_array, ('-G', '-C', -78))
        modify_single_weight(eps_array, ('G-', 'C-', -107))
    else:
        modify_single_weight(eps_array, ('AA', 'UU', -31))
        modify_single_weight(eps_array, ('AU', 'UA', -14))
        modify_single_weight(eps_array, ('AC', 'UG', -7))
        modify_single_weight(eps_array, ('AG', 'UC', -17))
        modify_single_weight(eps_array, ('UA', 'AU', -28))
        modify_single_weight(eps_array, ('UU', 'AA', -34))
        modify_single_weight(eps_array, ('UC', 'AG', -30))
        modify_single_weight(eps_array, ('UG', 'AC', -17))
        modify_single_weight(eps_array, ('CA', 'GU', -13))
        modify_single_weight(eps_array, ('CU', 'GA', -28))
        modify_single_weight(eps_array, ('CC', 'GG', -18))
        modify_single_weight(eps_array, ('CG', 'GC', -41))
        modify_single_weight(eps_array, ('GA', 'CU', -3))
        modify_single_weight(eps_array, ('GU', 'CA', -33))
        modify_single_weight(eps_array, ('GC', 'CG', -61))
        modify_single_weight(eps_array, ('GG', 'CC', -0))

        modify_single_weight(eps_array, ('-A', '-U', 106))
        modify_single_weight(eps_array, ('A-', 'U-', 102))
        modify_single_weight(eps_array, ('-U', '-A', 108))
        modify_single_weight(eps_array, ('U-', 'A-', 104))
        modify_single_weight(eps_array, ('-C', '-G', 121))
        modify_single_weight(eps_array, ('C-', 'G-', 98 ))
        modify_single_weight(eps_array, ('-G', '-C', 121))
        modify_single_weight(eps_array, ('G-', 'C-', 98 ))

def experimental_weights_37_deg(eps_array, neg=False, ):
    if not neg:
        modify_single_weight(eps_array, ('AA', 'UU', 40))
        modify_single_weight(eps_array, ('AU', 'UA', 14))
        modify_single_weight(eps_array, ('AC', 'UG', 1 ))
        modify_single_weight(eps_array, ('AG', 'UC', 29))
        modify_single_weight(eps_array, ('UA', 'AU', 37))
        modify_single_weight(eps_array, ('UU', 'AA', 43))
        modify_single_weight(eps_array, ('UC', 'AG', 26))
        modify_single_weight(eps_array, ('UG', 'AC', 17))
        modify_single_weight(eps_array, ('CA', 'GU', 9 ))
        modify_single_weight(eps_array, ('CU', 'GA', 34))
        modify_single_weight(eps_array, ('CC', 'GG', 18))
        modify_single_weight(eps_array, ('CG', 'GC', 48))
        modify_single_weight(eps_array, ('GA', 'CU', -1))
        modify_single_weight(eps_array, ('GU', 'CA', 39))
        modify_single_weight(eps_array, ('GC', 'CG', 84))
        modify_single_weight(eps_array, ('GG', 'CC', -8))

        modify_single_weight(eps_array, ('-A', '-U', -132))
        modify_single_weight(eps_array, ('A-', 'U-', -119))
        modify_single_weight(eps_array, ('-U', '-A', -124))
        modify_single_weight(eps_array, ('U-', 'A-', -126))
        modify_single_weight(eps_array, ('-C', '-G', -161))
        modify_single_weight(eps_array, ('C-', 'G-', -124))
        modify_single_weight(eps_array, ('-G', '-C', -81 ))
        modify_single_weight(eps_array, ('G-', 'C-', -130))
    else:
        modify_single_weight(eps_array, ('AA', 'UU', -40))
        modify_single_weight(eps_array, ('AU', 'UA', -14))
        modify_single_weight(eps_array, ('AC', 'UG', -1 ))
        modify_single_weight(eps_array, ('AG', 'UC', -29))
        modify_single_weight(eps_array, ('UA', 'AU', -37))
        modify_single_weight(eps_array, ('UU', 'AA', -43))
        modify_single_weight(eps_array, ('UC', 'AG', -26))
        modify_single_weight(eps_array, ('UG', 'AC', -17))
        modify_single_weight(eps_array, ('CA', 'GU', -9 ))
        modify_single_weight(eps_array, ('CU', 'GA', -34))
        modify_single_weight(eps_array, ('CC', 'GG', -18))
        modify_single_weight(eps_array, ('CG', 'GC', -48))
        modify_single_weight(eps_array, ('GA', 'CU', 1))
        modify_single_weight(eps_array, ('GU', 'CA', -39))
        modify_single_weight(eps_array, ('GC', 'CG', -84))
        modify_single_weight(eps_array, ('GG', 'CC', 8))

        modify_single_weight(eps_array, ('-A', '-U', 132))
        modify_single_weight(eps_array, ('A-', 'U-', 119))
        modify_single_weight(eps_array, ('-U', '-A', 124))
        modify_single_weight(eps_array, ('U-', 'A-', 126))
        modify_single_weight(eps_array, ('-C', '-G', 161))
        modify_single_weight(eps_array, ('C-', 'G-', 124))
        modify_single_weight(eps_array, ('-G', '-C', 81 ))
        modify_single_weight(eps_array, ('G-', 'C-', 130))

def experimental_ps_dna_rna(base, dsm_base, temp=50, transpose=True):
    eps_array = np.zeros((6, 6, 6, 6), dtype=np.int64) - base
    if temp == 50:
        experimental_weights_50_deg(eps_array, neg=False)
    elif temp == 37:
        experimental_weights_37_deg(eps_array, neg=False)
    if transpose:
        mini_transpose(eps_array)

    dsm_variable_to_csv(dsm_name=dsm_base, func=rna_dna_to_ps_dna_rna, kwargs={'eps_array': -eps_array},
                        func_extend=rna_dna_to_dna_rna, kwargs_extend=None)

def experimental_neg_ps_dna_rna(base, dsm_base, temp, transpose=True):
    eps_array = np.zeros((6, 6, 6, 6), dtype=np.int64) - base
    if temp == 50:
        experimental_weights_50_deg(eps_array, neg=True)
    elif temp == 37:
        experimental_weights_37_deg(eps_array, neg=True)
    if transpose:
        mini_transpose(eps_array)

    modify_single_weight(eps_array, ('-A', '-U', -106))
    modify_single_weight(eps_array, ('A-', 'U-', -102))
    modify_single_weight(eps_array, ('-U', '-A', -108))
    modify_single_weight(eps_array, ('U-', 'A-', -104))
    modify_single_weight(eps_array, ('-C', '-G', -121))
    modify_single_weight(eps_array, ('C-', 'G-', -98))
    modify_single_weight(eps_array, ('-G', '-C', -121))
    modify_single_weight(eps_array, ('G-', 'C-', -98))

    dsm_variable_to_csv(dsm_name=dsm_base, func=rna_dna_to_ps_dna_rna, kwargs={'eps_array': -eps_array},
                        func_extend=rna_dna_to_dna_rna, kwargs_extend=None)


def md_ps_dna_rna(base, dsm_base, punish_base=0):
    eps_array = np.zeros((6, 6, 6, 6), dtype=np.int64) - base
    modify_single_weight(eps_array, ('AA', 'UU', 114))
    modify_single_weight(eps_array, ('AU', 'UA', 103))
    modify_single_weight(eps_array, ('AC', 'UG', 17))
    modify_single_weight(eps_array, ('AG', 'UC', 174))
    modify_single_weight(eps_array, ('UA', 'AU', 132))
    modify_single_weight(eps_array, ('UU', 'AA', 149))
    modify_single_weight(eps_array, ('UC', 'AG', 191))
    modify_single_weight(eps_array, ('UG', 'AC', 160))
    modify_single_weight(eps_array, ('CA', 'GU', 177))
    modify_single_weight(eps_array, ('CU', 'GA', 85))
    modify_single_weight(eps_array, ('CC', 'GG', 207))
    modify_single_weight(eps_array, ('CG', 'GC', 215))
    modify_single_weight(eps_array, ('GA', 'CU', 199))
    modify_single_weight(eps_array, ('GU', 'CA', 151))
    modify_single_weight(eps_array, ('GC', 'CG', 218))
    modify_single_weight(eps_array, ('GG', 'CC', 271))

    modify_single_weight(eps_array, ('-A', '-U', -106 - punish_base))
    modify_single_weight(eps_array, ('A-', 'U-', -102 - punish_base))
    modify_single_weight(eps_array, ('-U', '-A', -108 - punish_base))
    modify_single_weight(eps_array, ('U-', 'A-', -104 - punish_base))
    modify_single_weight(eps_array, ('-C', '-G', -121 - punish_base))
    modify_single_weight(eps_array, ('C-', 'G-', -98 - punish_base))
    modify_single_weight(eps_array, ('-G', '-C', -121 - punish_base))
    modify_single_weight(eps_array, ('G-', 'C-', -98 - punish_base))

    dsm_variable_to_csv(dsm_name=dsm_base, func=rna_dna_to_ps_dna_rna, kwargs={'eps_array': -eps_array},
                        func_extend=rna_dna_to_dna_rna, kwargs_extend=None)
