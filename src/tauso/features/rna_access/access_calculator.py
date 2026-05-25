import numpy as np
import pandas as pd


def get_sense_with_flanks(pre_mrna: str, sense_start: int, sense_length: int, flank_size: int) -> str:
    """Sense window plus `flank_size` nt on each side, clipped to the sequence ends."""
    start = max(0, sense_start - flank_size)
    end = min(len(pre_mrna), sense_start + sense_length + flank_size)
    return pre_mrna[start:end]


def window_access_energies(access_res, rna_size, access_size, seed_sizes):
    """Per-position accessibility energy of an `access_size`-wide window, one column per seed.

    `access_res` is raccess' per-position opening energy for each seed size (one
    column per seed). For every window start position we sample seed-sized
    sub-windows spaced half a seed apart (step = seed // 2), take their weighted
    mean -- the final sub-window is clipped to the window edge and down-weighted by
    the fractional step -- and rescale by `access_size / seed`. raccess reports
    energies from the 3' end, hence the reversed lookup index.

    Returns a DataFrame with one `{seed}_avg` column per seed, indexed by window
    start position (0 .. rna_size - access_size).
    """
    positions = np.arange(rna_size - access_size + 1)

    avg_cols = {}
    for seed in seed_sizes:
        step = seed // 2
        offsets = np.append(np.arange(0, access_size - seed, step), access_size - seed)

        weights = np.ones(len(offsets))
        if len(offsets) > 1:
            weights[-1] = (offsets[-1] - offsets[-2]) / step

        energies = access_res[seed].to_numpy()
        abs_offsets = offsets[None, :] + positions[:, None]  # (n_pos, n_offsets)
        sampled = energies[len(energies) - 1 - abs_offsets]  # 3'->5' lookup

        scaled = sampled * (access_size / seed) * weights
        avg_cols[f"{seed}_avg"] = scaled.sum(axis=1) / weights.sum()

    return pd.DataFrame(avg_cols, index=positions)
