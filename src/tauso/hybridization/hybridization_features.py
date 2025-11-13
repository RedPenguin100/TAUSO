from typing import List

# https://pmc.ncbi.nlm.nih.gov/articles/PMC9116672/#ack1

# DNA/DNA - PSDNA/DNA at 50 degrees
exp_ps_diff_weights_50 = {
    'AA': -31,
    'AU': -14,
    'AC': -6,
    'AG': -23,
    'UA': -28,
    'UU': -34,
    'UC': -24,
    'UG': -18,
    'CA': -13,
    'CU': -28,
    'CC': -18,
    'CG': -41,
    'GA': -3,
    'GU': -33,
    'GC': -61,
    'GG': -0,
}

# DNA/DNA - PSDNA/DNA at 37 degrees
# calculated using (enthlapy_mean * 1000 - 310 * entropy_mean) / 1000

exp_ps_diff_weights_37 = {
    'AA': -40,
    'AU': -14,
    'AC': -1,
    'AG': -30,
    'UA': -37,
    'UU': -34,
    'UC': -26,
    'UG': -17,
    'CA': -10,
    'CU': -34,
    'CC': -18,
    'CG': -48,
    'GA': 2,
    'GU': -39,
    'GC': -84,
    'GG': 8,
}

# DNA/RNA weights
# https://pubs.acs.org/doi/10.1021/bi00035a029
exp_rna_weights_37 = {
    'UU': -100,
    'GU': -210,
    'CU': -180,
    'AU': -90,
    'UG': -90,
    'GG' : -210,
    'CG' : -170,
    'AG' : -90,
    'UC' : -130,
    'GC' : -270,
    'CC' : -290,
    'AC' : -110,
    'UA' : -60,
    'GA' : -150,
    'CA' : -160,
    'AA' : -20
}


def get_exp_psrna_hybridization(seq: str, temp=37) -> float:
    total_hybridization = 0
    for i in range(len(seq) - 1):
        L, R = seq[i], seq[i + 1]
        # if temp == 50:
        #     total_hybridization += exp_rna_weights_37[L+R] - exp_ps_diff_weights_50[L + R]
        if temp == 37:
            total_hybridization += exp_rna_weights_37[L+R] - exp_ps_diff_weights_37[L + R]
        else:
            total_hybridization = 0
    return total_hybridization


def get_exp_psrna_hybridization_normalized(seq: str, temp=50) -> float:
    return get_exp_psrna_hybridization(seq, temp) / len(seq)
