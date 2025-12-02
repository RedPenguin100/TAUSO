import pytest

from tauso.genome.read_human_genome import get_locus_to_data_dict


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