def calc_tm_14nt(bases_counter):
    return 64.9 + 41.0 * ((bases_counter['G'] + bases_counter['C'] - 16.4) / (
                bases_counter['G'] + bases_counter['C'] + bases_counter['A'] + bases_counter['T']))


def calc_tm_13nt(bases_counter):
    return (bases_counter['A'] + bases_counter['T']) * 2 + (bases_counter['G'] + bases_counter['C']) * 4


def calc_tm(seq):
    seq_upper = seq.upper()
    bases_counter = {nt: seq_upper.count(nt) for nt in "ATGC"}
    if len(seq) > 13:
        return calc_tm_14nt(bases_counter)
    else:
        return calc_tm_13nt(bases_counter)
