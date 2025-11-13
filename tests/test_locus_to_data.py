import pytest
import pickle

from tauso.consts import CACHE_DIR
from tauso.timer import Timer
from tauso.read_yeast import get_locus_to_data_dict_alternative, get_locus_to_data_dict_yeast
from tauso.read_human_genome import get_locus_to_data_dict


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


@pytest.mark.slow
def test_intron_regression():
    intron_yeast_gene = 'YNCA0002W'

    locus_to_data = get_locus_to_data_dict_alternative()

    assert len(locus_to_data[intron_yeast_gene].introns) == 1
    intron = locus_to_data[intron_yeast_gene].introns[0]
    assert str(intron) == 'CGACTTCCTGATTAAACAGGAAGACAAAGCA'


@pytest.mark.slow
def test_begin_regression_human():
    locus_to_data = get_locus_to_data_dict(gene_subset=['KLKB1'])

    first_exon_begin = locus_to_data['KLKB1'].exons[0][0:33]

    assert str(first_exon_begin) == 'AGTGCCACATTAGAACAGCTTGAAGACCGTTCA'

    full_mrna_begin = locus_to_data['KLKB1'].full_mrna
    idx = full_mrna_begin.find(first_exon_begin)
    assert idx == 18528

# not slow!
def test_begin_regression_fast_human():
    cache_path = CACHE_DIR / 'gene_to_data_simple_cache.pickle'
    with open(cache_path, 'rb') as f:
        gene_to_data = pickle.load(f)

    first_exon_begin = gene_to_data['KLKB1'].exons[0][0:33]
    assert str(first_exon_begin) == 'AGTGCCACATTAGAACAGCTTGAAGACCGTTCA'
    full_mrna_begin = gene_to_data['KLKB1'].full_mrna
    idx = full_mrna_begin.find(first_exon_begin)
    assert idx == 18528
