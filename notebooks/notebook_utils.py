import pickle

import numpy as np

from notebooks.consts import INHIBITION, CACHE_DIR
from tauso.genome.read_human_genome import get_locus_to_data_dict


def log_correction(df, correction=0.01):
    df.loc[:, 'log_inhibition'] = -np.log(-df[INHIBITION] + 100 + correction)  # to avoid log 0


def read_cached_gene_to_data(genes_u):
    # Load gene information from cache or generate it if needed
    cache_path = CACHE_DIR / 'gene_to_data_simple_cache.pickle'
    if not cache_path.exists():
        gene_to_data = get_locus_to_data_dict(include_introns=True, gene_subset=genes_u)
        with open(cache_path, 'wb') as f:
            pickle.dump(gene_to_data, f)
    else:
        with open(cache_path, 'rb') as f:
            gene_to_data = pickle.load(f)
    return gene_to_data
