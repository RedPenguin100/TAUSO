import pandas as pd

from Bio import SeqIO
from fuzzysearch import find_near_matches
from numba import njit

from tauso.util import get_antisense
from tests.common.consts import INTEGRATION_TESTS_PATH


@njit
def iterate_template(template_seq, l_values):
    print("Hello")
    for l in l_values:
        for i in range(len(template_seq) - l + 1):
            yield i, l, template_seq[i:i + l]


@njit
def iterate_template_antisense(template_seq, l_values):
    print("Anti hello")
    for l in l_values:
        for i in range(len(template_seq) - l + 1):
            yield i, l, get_antisense(template_seq[i:i + l])


def get_gfp_seq_and_context():
    gfp_context_path = INTEGRATION_TESTS_PATH / 'data' / 'GFP_context.txt'
    gfp_first_exp_path = INTEGRATION_TESTS_PATH / 'data' / 'GFP_first_exp.fasta'

    gfp_obj = next(SeqIO.parse(str(gfp_first_exp_path), 'fasta'))
    gfp_seq = str(gfp_obj.seq.upper())

    with open(str(gfp_context_path), 'r') as f:
        gfp_context = f.read().upper()

    gfp_start = gfp_context.find(gfp_seq)
    if gfp_start == -1:
        raise ValueError("Context not found!")

    return (gfp_seq, gfp_context)


def get_gfp_first_exp(gap=100):
    # TODO: gap should be always 100 in this function
    gfp_seq, gfp_context = get_gfp_seq_and_context()

    gfp_start = gfp_context.find(gfp_seq)
    if gfp_start == -1:
        raise ValueError("Context not found!")

    gfp_ext = gfp_context[gfp_start - gap: gfp_start + len(gfp_seq) + gap]

    return gfp_ext


def get_gfp_second_exp():
    right_gap = 50
    gfp_seq, gfp_context = get_gfp_seq_and_context()

    gfp_start = gfp_context.find(gfp_seq)
    gfp_ext = gfp_context[gfp_start: gfp_start + len(gfp_seq) + right_gap]

    return gfp_ext


DEGRON_LENGTH = 126


def get_degron_and_gap_third_exp():
    left_gap = 25

    # Very specific for this sequence protein
    gfp_seq, gfp_context = get_gfp_seq_and_context()
    gfp_start = gfp_context.find(gfp_seq)
    degron_first = gfp_start + len(gfp_seq)  # one after gfp end
    degron_and_gap = gfp_context[degron_first - left_gap:degron_first + DEGRON_LENGTH]
    return degron_and_gap


# Includes - GFP, Degron, WPRE, 3UTR
def get_extended_gfp():
    gfp_seq, gfp_context = get_gfp_seq_and_context()
    gfp_start = gfp_context.find(gfp_seq)
    return gfp_context[gfp_start:]


def get_degron_gfp_scrambled_third_exp():
    gfp_seq, gfp_context = get_gfp_seq_and_context()
    gfp_start = gfp_context.find(gfp_seq)
    degron_end = gfp_start + len(gfp_seq) + DEGRON_LENGTH
    return gfp_context[gfp_start:degron_end]


def get_3utr_gfp():
    gfp_seq, gfp_context = get_gfp_seq_and_context()
    utr_start = 2813
    utr_length = 234

    three_utr = gfp_context[utr_start:utr_start + utr_length]
    return three_utr


def get_angtpl2():
    with open(DATA_PATH / 'ANGPTL2_pre_mrna.txt', 'r') as file:
        return file.read()

def get_bcl2_patent():
    with open(DATA_PATH / 'BCL2_human_premrna.txt', 'r') as file:
        return file.read()


def generate_scrambled(target_seq):
    l_values = [17, 18, 19, 20, 21]

    matches_per_distance = [0, 0, 0, 0]
    df = pd.DataFrame(
        columns=['sense_start', 'sense_length', '0_matches', '1_matches', '2_matches', '3_matches'])

    for i, l, antisense in iterate_template_antisense(target_seq, l_values):
        matches = find_near_matches(antisense, target_seq, max_insertions=0, max_deletions=0, max_l_dist=3)
        for match in matches:
            matches_per_distance[match.dist] += 1

        if (i, l) not in df.index:
            df.loc[len(df)] = [i, l, matches_per_distance[0], matches_per_distance[1], matches_per_distance[2],
                               matches_per_distance[3]]

    return df


def generate_scrambled_permutation(target_seq: str, l_values=[16,17,18,19,20,21]):
    parsing_type = '2'

    aso_to_scores = dict()

    for (i, l, sense) in iterate_template(target_seq, l_values):
        tmp_results = get_trigger_mfe_scores_by_risearch(sense, {'GFP': target_seq}, minimum_score=900,
                                                         neighborhood=l,
                                                         parsing_type=parsing_type)
        scores = get_mfe_scores(tmp_results, parsing_type)

        aso_to_scores[(i, l)] = sum(scores[0])

    return aso_to_scores
