import pickle

import numpy as np

from notebooks.consts import CACHE_DIR
from tauso.genome.read_human_genome import get_locus_to_data_dict
from tauso.data.consts import INHIBITION


def log_correction(df, correction=0.01):
    df.loc[:, 'log_inhibition'] = -np.log(-df[INHIBITION] + 100 + correction)  # to avoid log 0
