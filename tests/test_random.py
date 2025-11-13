import pytest
import ViennaRNA as RNA

from tauso.random_util import generate_random_dna
from tauso.features.seq_features import get_gc_content


@pytest.mark.parametrize('length', [16, 17, 18, 19, 20, 21, 22])
def test_random_dna(length):
    gc_lower = 0.5
    gc_upper = 0.65
    min_fold_energy = -1
    attempts = 100
    random_dna = generate_random_dna(length, gc_lower, gc_upper, min_fold_energy, attempts)

    print(len(random_dna))
    print(len(set(random_dna)))

    for dna in random_dna:
        gc_content = get_gc_content(dna)
        assert gc_content >= gc_lower
        assert gc_content <= gc_upper
        assert min_fold_energy <= RNA.fold(dna)[1]
        assert length == len(dna)
