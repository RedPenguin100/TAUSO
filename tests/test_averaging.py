import pytest

import pandas as pd
import numpy as np

from tauso.features.fold import get_weighted_energy


def test_window_average():
    energies = np.array([0, 1, 2, 3, 4, 5])

    step_size = 15
    window_size = 40

    l_values = [1]
    target_seq = np.zeros(step_size * len(energies))

    results = []
    columns = ['sense_start', 'sense_length', 'mean_fold']

    for l in l_values:
        for i in range(len(target_seq) - l + 1):
            target_start = i
            mean_fold = get_weighted_energy(target_start, l, step_size, energies, window_size)
            results.append((i, l, mean_fold))

    df = pd.DataFrame(results, columns=columns)

    assert df.loc[0, 'mean_fold'] == 0
    assert df.loc[14, 'mean_fold'] == 0
    assert df.loc[15, 'mean_fold'] == 0.5
    assert df.loc[29, 'mean_fold'] == 0.5
    assert df.loc[30, 'mean_fold'] == 1.
    assert df.loc[39, 'mean_fold'] == 1.
    assert df.loc[40, 'mean_fold'] == 1.5
    assert df.loc[44, 'mean_fold'] == 1.5
    assert df.loc[45, 'mean_fold'] == 2


