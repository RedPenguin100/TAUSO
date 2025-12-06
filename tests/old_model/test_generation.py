import pytest

from tauso.old_model_generation.main import get_target_sequence, get_init_df


@pytest.mark.integration
def test_stub_gene():
    exogenous_gene = 'GAUACAGAUACAGAUACA'
    assert exogenous_gene == get_target_sequence('AAA', exogenous_gene)


@pytest.mark.integration
def test_generate_candidates():
    short_gene_name = 'DDX11L1'
    short_gene =  get_target_sequence(short_gene_name)
    candidate_asos = get_init_df(short_gene, 'moe')
    print("Length: ", len(candidate_asos))
    print(candidate_asos)

