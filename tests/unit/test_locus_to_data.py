import pytest
import pickle

from tauso.consts import CACHE_DIR
from tauso.timer import Timer
from tauso.read_yeast import get_locus_to_data_dict_alternative, get_locus_to_data_dict_yeast
from tauso.genome.read_human_genome import get_locus_to_data_dict


# TODO: should CI run this?
@pytest.mark.slow
def test_locus():
    with Timer() as t:
        locus_to_data = get_locus_to_data_dict_yeast()
    print(f"Regular Took: {t.elapsed_time}s")

    with Timer() as t:
        locus_to_data_alt = get_locus_to_data_dict_alternative()
    print(f"DB Took: {t.elapsed_time}s")

    for key, value in locus_to_data.items():
        alt_exons = locus_to_data_alt[key].exons
        exons = locus_to_data[key].exons
        assert len(exons) == len(alt_exons)

        for i in range(len(exons)):
            exon = exons[i]
            alt_exon = alt_exons[i]

            assert exon == alt_exon, "key: " + key + " len exons " + str(len(exons)) + "i: " + str(i)


# TODO: decide what to do with yeast
# def test_intron_regression():
#     intron_yeast_gene = 'YNCA0002W'
#
#     locus_to_data = get_locus_to_data_dict_alternative()
#
#     assert len(locus_to_data[intron_yeast_gene].introns) == 1
#     intron = locus_to_data[intron_yeast_gene].introns[0]
#     assert str(intron) == 'CGACTTCCTGATTAAACAGGAAGACAAAGCA'


# TODO: this is not slow, but CI can't run this.
@pytest.mark.slow
def test_begin_exon_regression_human():
    locus_to_data = get_locus_to_data_dict(gene_subset=['KLKB1'])

    full_mrna_begin = locus_to_data['KLKB1'].full_mrna
    idx = full_mrna_begin.find('AGTGCCACATTAGAACAGCTTGAAGACCGTTCA')
    print("Idx actual: ", idx)
    assert idx == 2083


# def test_begin_intron_regression_human():
#     locus_to_data = get_locus_to_data_dict(gene_subset=['KLKB1'])
#
#     first_exon_begin = locus_to_data['KLKB1'].exons[0][0:33]
#
#     assert str(first_exon_begin) == 'AGTGCCACATTAGAACAGCTTGAAGACCGTTCA'
#
#     full_mrna_begin = locus_to_data['KLKB1'].full_mrna
#     idx = full_mrna_begin.find(first_exon_begin)
#     assert idx == 18528