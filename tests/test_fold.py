import pytest
from Bio import SeqIO

from tauso.consts import DATA_PATH
# from external.risearch.RIsearch1.numpy_to_csv import dsm_variable_to_csv,numpy_to_csv

from tauso.features.vienna_fold import get_mfe_scores
from tauso.hybridization.fast_hybridization import get_trigger_mfe_scores_by_risearch, Interaction, dump_target_file
from tauso.util import get_antisense


def get_gfp_seq_and_context():
    gfp_context_path = DATA_PATH / 'GFP_context.txt'
    gfp_first_exp_path = DATA_PATH / 'GFP_first_exp.fasta'

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


@pytest.mark.skip
def test_empty():
    name_to_sequence = {"GFP": get_gfp_first_exp()}

    result = get_trigger_mfe_scores_by_risearch(name_to_sequence['GFP'][0:20], name_to_sequence, minimum_score=30000,
                                                # very high score, no candidates will be found
                                                parsing_type='2')
    mfe_scores = get_mfe_scores(result, '2')

    assert len(mfe_scores) == 1
    assert len(mfe_scores[0]) == 0


@pytest.mark.skip
def test_risearch():
    # example
    name_to_sequence = {
        "0": "AUGUGUUCAAUUUUUGGCGUAUUCGAUAUCAAAACAGACGCAGUUGAGCUGCGUAAGAAAGCACUCGAGCUGUCACGCCUGAUGCGUCAUCGUGGCCCGGACUGGUCCGGUAUUUAUGCCAGCGAUAACGCCAUUCUCGCCCACGAACGUCUGUCAAUUGUUGACGUUAACGCGGGGGCGCAACCUCUCUACAACCAACAAAAAACUCAUGUGCUGGCGGUAAACGGUGAAAUCUACAACCACCAGGCACUGCGCGCCGAAUAUGGCGAUCGUUAUCAGUUCCAGACCGGAUCUGACUGUGAAGUGAUCCUCGCGCUGUAUCAGGAAAAAGGGCCGGAAUUUCUCGACGACUUGCAGGGCAUGUUUGCCUUUGCCUUGUACGACAGCGAAAAAGAUGCCUACCUGAUUGGUCGCGACCAUCUGGGGAUCAUCCCACUGUAUAUGGGCUAUGACGAACACGGUCAGCUGUAUGUGGCCUCAGAAAUGAAAGCCCUGGUGCCAGUUUGCCGCACGAUUAAAGAGUUCCCGGCGGGGAGCUAUUUGUGGAGCCAGGACGGCGAAAUCCGUUCUUAUUAUCAUCGCGACUGGUUCGACUACGAUGCGGUGAAAGAUAACGUAACCGACAAAAACGAGCUGCGUCAGGCACUUGAAGAUUCCGUUAAAAGCCAUCUGAUGUCUGAUGUGCCUUACGGUGUGCUGCUUUCUGGUGGUCUGGAUUCCUCAAUUAUUUCCGCUAUCACCAAGAAAUACGCAGCCCGUCGCGUGGAAGAUCAGGAACGCUCUGAAGCCUGGUGGCCGCAGUUACACUCCUUUGCUGUAGGUCUGCCGGGUUCACCGGAUCUUAAAGCAGCCCAGGAAGUGGCAAACCAUCUGGGCACGGUGCAUCACGAAAUUCACUUCACUGUACAGGAAGGUCUGGAUGCCAUCCGCGACGUGAUUUACCACAUCGAAACUUAUGAUGUGACCACAAUUCGCGCUUCAACACCGAUGUAUUUAAUGUCGCGUAAGAUCAAGGCGAUGGGCAUUAAAAUGGUGCUGUCCGGUGAAGGUUCUGAUGAAGUGUUUGGCGGUUAUCUUUACUUCCAUAAAGCGCCCAACGCUAAAGAACUGCAUGAAGAGACGGUGCGUAAACUGCUGGCCCUGCAUAUGUAUGACUGCGCGCGCGCCAACAAAGCGAUGUCAGCCUGGGGCGUGGAAGCACGCGUUCCGUUCCUCGACAAAAAAUUCCUUGAUGUGGCGAUGCGCAUUAACCCGCAGGAUAAAAUGUGCGGUAACGGCAAAAUGGAAAAACACAUCCUGCGUGAAUGUUUUGAGUCAUACCUGCCCGCAAGCGUGGCCUGGCGGCAGAAAGAGCAGUUCUCCGAUGGCGUCGGUUACAGUUGGAUCGACACCCUGAAAGAAGUGGCGGCGCAGCAGGUUUCUGAUCAGCAACUGGAAACUGCCCGCUUCCGCUUCCCGUACAACACGCCAACCUCAAAAGAAGCGUAUCUGUACCGGGAGAUCUUUGAAGAACUGUUCCCGCUUCCGAGCGCCGCUGAGUGCGUGCCUGGCGGUCCUUCCGUCGCGUGUUCUUCCGCUAAAGCGAUCGAAUGGGAUGAAGCGUUCAAGAAAAUGGACGAUCCAUCUGGUCGUGCGGUUGGUGUUCACCAGUCGGCAUAUAAGUAA",
        "1": "AUGUUCGAACAACGCGUAAAUUCUGACGUACUGACCGUUUCUACCGUUAACUCUCAGGAUCAGGUAACCCAAAAACCCCUGCGUGACUCGGUUAAACAGGCACUGAAGAACUAUUUUGCUCAACUGAAUGGUCAGGAUGUGAAUGACCUCUAUGAGCUGGUACUGGCUGAAGUAGAACAGCCCCUGUUGGACAUGGUGAUGCAAUACACCCGUGGUAACCAGACCCGUGCUGCGCUGAUGAUGGGCAUCAACCGUGGUACGCUGCGUAAAAAAUUGAAAAAAUACGGCAUGAACUAA"
    }

    target_path = dump_target_file('target-cache.fa', name_to_sequence)
    result = get_trigger_mfe_scores_by_risearch("UAGAUGCGCCACUUGUGGUAUUCCCGCAUC", name_to_sequence, minimum_score=900,
                                                parsing_type='2', target_file_cache=str(target_path))
    print(result)
    mfe_scores = get_mfe_scores(result, '2')

    print(mfe_scores)

@pytest.mark.skip
def test_bad_fit():
    sense = 'TTTTTTTCTTCCATT'

    result = get_trigger_mfe_scores_by_risearch(sense, {'0': sense + sense}, parsing_type='2', minimum_score=900)
    print(result)

    gfp_seq = get_gfp_first_exp(gap=0)
    sample_seq = gfp_seq[0:15]
    result = get_trigger_mfe_scores_by_risearch('T' + sample_seq[1:15], {"gfp": gfp_seq}, parsing_type='2',
                                                minimum_score=900)
    print(result)


@pytest.mark.skip
def test_risearch_gfp():
    gfp_seq = get_gfp_first_exp(gap=0)
    sample_seq = gfp_seq[:20]
    print("GFP ontarget(?) : ", gfp_seq[695:714])
    print("Sample", sample_seq)
    print("Sample antisense", get_antisense(sample_seq))
    # name_to_seq = {f"gfp_seq{i}" : gfp_seq for i in range(100)}
    name_to_seq = {f"gfp_seq": gfp_seq}
    result = get_trigger_mfe_scores_by_risearch(sample_seq, name_to_seq,
                                                interaction_type=Interaction.DNA_RNA_NO_WOBBLE,
                                                minimum_score=900, neighborhood=30, parsing_type='2')
    print(result)

    mfe_scores = get_mfe_scores(result, '2')
    print(mfe_scores)

    bad_samples = [s + sample_seq[3:20] for s in ["AAA", "ATA", "AGA", "ACG"]]

    for bad_sample in bad_samples:
        result = get_trigger_mfe_scores_by_risearch(bad_sample, name_to_seq,
                                                    interaction_type=Interaction.DNA_RNA_NO_WOBBLE, minimum_score=900,
                                                    neighborhood=30, parsing_type='2')
        print(result)

        mfe_scores = get_mfe_scores(result, '2')
        print(mfe_scores)


# from external.risearch.RIsearch1.numpy_to_csv import dsm_variable_to_csv,numpy_to_csv
# from tauso.modified_dsm import make_dsm_ps_dna_rna , make_dsm_dna_rna

# TODO: uncomment when RIsearch is integrated
# def test_risearch_gfp_modified():
#     dsm_variable_to_csv()
#     make_dsm_dna_rna()
#
#     gfp_seq = get_gfp_first_exp(gap=0)
#     sample_seq = gfp_seq[:20]
#     print("GFP ontarget(?) : ", gfp_seq[695:714])
#     print("Sample", sample_seq)
#     print("Sample antisense", get_antisense(sample_seq))
#     # name_to_seq = {f"gfp_seq{i}" : gfp_seq for i in range(100)}
#     name_to_seq = {f"gfp_seq": gfp_seq}
#     result = get_trigger_mfe_scores_by_risearch(sample_seq, name_to_seq,
#                                                 interaction_type=Interaction.MODIFIED,
#                                                 minimum_score=900, neighborhood=30, parsing_type='2')
#     print(result)
#
#     mfe_scores = get_mfe_scores(result, '2')
#     print(mfe_scores)
#
#     bad_samples = [s + sample_seq[3:20] for s in ["AAA", "ATA", "AGA", "ACG"]]
#
#     for bad_sample in bad_samples:
#         result = get_trigger_mfe_scores_by_risearch(bad_sample, name_to_seq,
#                                                     interaction_type=Interaction.MODIFIED, minimum_score=900,
#                                                     neighborhood=30, parsing_type='2')
#         print(result)
#
#         mfe_scores = get_mfe_scores(result, '2')
#         print(mfe_scores)