import math
import random
import ViennaRNA as RNA


def quasi_normal_random_int(a: int, b: int) -> int:
    mean = (a + b) / 2
    std_dev = (b - a) / 6

    normal_num = random.gauss(mean, std_dev)
    if a <= normal_num <= b:
        return int(round(normal_num))
    if normal_num > b:
        return b
    return a


def generate_random_dna(length, gc_lower=0.5, gc_upper=0.65, min_fold_energy=-1., attempts=100):
    seqs = []
    for attempt in range(attempts):
        gc_lower_amount = math.ceil(gc_lower * length)
        gc_upper_amount = math.floor(gc_upper * length)
        if gc_lower_amount > gc_upper_amount:
            raise ValueError(f"gc_content {gc_lower} <= X <= {gc_upper_amount} doesn't exist for length {length}")

        # gc content generation is uniform within boundaries
        gc_generated_amount = random.randint(gc_lower_amount, gc_upper_amount)
        at_generated_amount = length - gc_generated_amount
        # Nucleotide generation is normal to minimize skewed A:T or G:C ratios
        g_amount = quasi_normal_random_int(0, gc_generated_amount)
        c_amount = gc_generated_amount - g_amount
        a_amount = quasi_normal_random_int(0, at_generated_amount)
        t_amount = at_generated_amount - a_amount

        base_dna = ["G"] * g_amount + ["C"] * c_amount + ["A"] * a_amount + ["T"] * t_amount
        random.shuffle(base_dna)
        dna_string = ''.join(base_dna)
        if RNA.fold(dna_string)[1] > min_fold_energy:
            seqs.append(dna_string)

    return list(set(seqs))
