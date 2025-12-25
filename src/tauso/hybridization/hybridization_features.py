from typing import List

from tauso.hybridization.weights.dna import DNA_DNA_WEIGHTS
from tauso.hybridization.weights.lna import LNA_DNA_WEIGHTS
from tauso.util import get_nucleotide_watson_crick

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
    'GG': -210,
    'CG': -170,
    'AG': -90,
    'UC': -130,
    'GC': -270,
    'CC': -290,
    'AC': -110,
    'UA': -60,
    'GA': -150,
    'CA': -160,
    'AA': -20
}


def get_exp_dna_rna_hybridization(seq: str, temp=37) -> float:
    seq = seq.replace('T', 'U')

    total_hybridization = 0
    for i in range(len(seq) - 1):
        L, R = seq[i], seq[i + 1]
        if temp == 37:
            total_hybridization += exp_rna_weights_37[L + R]
        else:
            total_hybridization = 0
    return total_hybridization


def get_exp_psrna_hybridization(seq: str, temp=37) -> float:
    seq = seq.replace('T', 'U')

    total_hybridization = 0
    for i in range(len(seq) - 1):
        L, R = seq[i], seq[i + 1]
        if temp == 37:
            total_hybridization += exp_rna_weights_37[L + R] - exp_ps_diff_weights_37[L + R]
        else:
            total_hybridization = 0
    return total_hybridization


def get_exp_psrna_hybridization_normalized(seq: str, temp=50) -> float:
    return get_exp_psrna_hybridization(seq, temp) / len(seq)


# --- Usage ---
# Replace 'ACTIVITY_SCORE' with your actual target column name (e.g., 'Kd', 'Tm', 'Inhibition')
# target_column = 'ACTIVITY_SCORE'
def calculate_3rd_gen_diff(seq_3to5, fmt_3to5, lna_params, temp_c=37.0, letter='L'):
    # Reverse to 5'->3' for processing
    seq = seq_3to5.upper()[::-1]
    fmt = fmt_3to5.upper()[::-1]

    if len(seq) != len(fmt): return None

    temp_k = temp_c + 273.15

    total_dH = 0.0
    total_dS = 0.0

    # --- 1. Add Initiation Penalty (Once per strand) ---
    # total_dH += INITIATION['dH']
    # total_dS += INITIATION['dS']

    # --- 2. Iterate Stacks ---
    for i in range(len(seq) - 1):
        b1, b2 = seq[i], seq[i + 1]
        m1, m2 = fmt[i], fmt[i + 1]

        # Determine Lookup Key
        if m1 == 'd' and m2 == 'd':
            continue
        else:
            # CASE B: LNA Involved (Use Owczarzy/McTigue)
            if m1 == letter and m2 == letter:
                top = f"+{b1}+{b2}"
            elif m1 == 'd' and m2 == letter:
                top = f"{b1}+{b2}"
            elif m1 == letter and m2 == 'd':
                top = f"+{b1}{b2}"
            else:
                continue

            c1, c2 = get_nucleotide_watson_crick(b1), get_nucleotide_watson_crick(b2)
            key = f"{top}/{c1}{c2}"

            if key in lna_params:
                total_dH += lna_params[key]['dH']
                total_dS += lna_params[key]['dS']
            else:
                print("Whoops")

    # Final dG calculation
    total_dG = total_dH - (temp_k * (total_dS / 1000.0))
    return round(total_dG, 3)


def calculate_lna(antisense, chemical_pattern):
    return calculate_3rd_gen_diff(
        antisense,
        chemical_pattern,
        LNA_DNA_WEIGHTS,  # Your existing LNA Dict
    )

def calculate_cet(antisense, chemical_pattern):
    return calculate_3rd_gen_diff(
        antisense,
        chemical_pattern,
        LNA_DNA_WEIGHTS,  # Your existing LNA Dict
        letter='C'
    )



def calculate_dna(antisense, temp_c=37.0):
    # Reverse to 5'->3' for processing
    seq = antisense.upper()[::-1]
    seq = seq.replace('U', 'T')

    temp_k = temp_c + 273.15

    total_dH = 0.0
    total_dS = 0.0

    # --- 1. Add Initiation Penalty (Once per strand) ---
    # total_dH += INITIATION['dH']
    # total_dS += INITIATION['dS']

    # --- 2. Iterate Stacks ---
    for i in range(len(seq) - 1):
        b1, b2 = seq[i], seq[i + 1]


        top = f"{b1}{b2}"
        c1, c2 = get_nucleotide_watson_crick(b1), get_nucleotide_watson_crick(b2)
        key = f"{top}/{c1}{c2}"

        if key in DNA_DNA_WEIGHTS:
            total_dH += DNA_DNA_WEIGHTS[key]['dH']
            total_dS += DNA_DNA_WEIGHTS[key]['dS']
        else:
            print("Whoops")

    # Final dG calculation
    total_dG = total_dH - (temp_k * (total_dS / 1000.0))
    return round(total_dG, 3)


def calc_methylcytosines(sequence, chemical_pattern, modification):
    if 'methylcytosines' not in modification:
        return 0
    gap_count = 0
    for base, char in zip(sequence, chemical_pattern):
        if char == 'd':
            if base == 'C':
                gap_count +=1
    return gap_count