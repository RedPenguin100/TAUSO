from typing import NamedTuple

import numpy as np

from ...util import dna_to_rna
from ._vienna_internal import window_mfe


class SharedFold(NamedTuple):
    """The one folded sequence a sub-sequence's windows are read out of.

    `sequence` is the sequence that was folded, `offset` is where the sub-sequence being
    swept starts inside it, and `max_bp_span` is the widest window that will be asked
    for, so several window sizes land on the same fold rather than one fold each.
    """

    sequence: str
    offset: int
    max_bp_span: int


def _get_cached_mfe(sequence, start=0, window_size=None, fold=None):
    """MFE of ``sequence[start:start + window_size]``, sharing one fold of `sequence`.

    The whole `sequence` is folded once and cached, and each window is read back out of
    those shared matrices rather than folded on its own. Pass a `fold` to read the
    window out of some wider sequence that `sequence` sits inside instead.
    """
    if window_size is None:
        window_size = len(sequence) - start
    if fold is None:
        fold = SharedFold(sequence, 0, window_size)
    return window_mfe(fold.sequence, fold.offset + start, window_size, fold.max_bp_span)


def _per_position_avg_energies(sequence, window_size, start, stop, step, fold=None):
    """One sliding-window sweep over `sequence`; returns the per-position average MFE.

    Fold each window in range(start, stop, step), spread its mfe/window_size
    evenly over the window's positions, then average per position. Positions no
    window covers are NaN.
    """
    seq_len = len(sequence)
    energy_values = np.zeros(seq_len)
    counts = np.zeros(seq_len)
    for i in range(start, stop, step):
        # get_catched_mfe is equivalent to RNA.fold(sequence[i:i+window_size])[1], but faster.
        # only fold the entire sequence once, and then fetch folded items.
        mfe_per_nt = _get_cached_mfe(sequence, i, window_size, fold) / window_size
        energy_values[i : i + window_size] += mfe_per_nt
        counts[i : i + window_size] += 1
    return np.divide(
        energy_values,
        counts,
        out=np.full_like(energy_values, np.nan),
        where=counts != 0,
    )


def calculate_avg_mfe_per_step(sequence, sense_start_in_flank, sense_length, window_size, steps, fold=None):
    """Sliding-window MFE averaged over the sense region, for several `step` values at once.

    The plain version of each step sweeps windows from position 0 across the whole
    sequence and nanmeans the sense-region positions. Two shortcuts give the same
    numbers more cheaply:
      1. Only sweep windows overlapping the sense region — outside windows touch
         only positions the final mean never reads.
      2. Overlapping windows share a single fold, each reading its own energy back
         out of the shared matrices.
    """
    sequence = dna_to_rna(sequence)
    seq_len = len(sequence)
    sense_end = sense_start_in_flank + sense_length

    if not (0 <= sense_start_in_flank < seq_len and sense_end <= seq_len):
        return {step: np.nan for step in steps}

    # A window at position i overlaps the sense region iff
    # i+window-1 >= sense_start AND i <= sense_end-1.
    min_i = max(0, sense_start_in_flank - window_size + 1)
    max_i = min(seq_len - window_size, sense_end - 1)

    out = {}
    for step in steps:
        # First multiple of `step` >= min_i. Same grid as starting at 0, just filtered.
        start_i = ((min_i + step - 1) // step) * step
        avg_energies = _per_position_avg_energies(sequence, window_size, start_i, max_i + 1, step, fold)
        out[step] = np.nanmean(avg_energies[sense_start_in_flank:sense_end])
    return out


def calculate_avg_mfe_per_setting(mrna, global_start, sense_length, settings, fold_region=None):
    """Run `calculate_avg_mfe_per_step` for every (flank, window, step) in `settings`.

    Each setting cuts its own flank-padded sub-sequence around the target site. Those
    cuts are all substrings of the widest one, so the widest is folded once and every
    setting reads its windows out of it. Returns a dict keyed by the setting triples.

    `fold_region` is the (start, end) slice of `mrna` to fold. It defaults to exactly the
    widest cut; pass a wider region shared with neighbouring target sites and they all
    read their windows out of that one fold.
    """
    widest_flank = max(flank for flank, _, _ in settings)
    widest_window = max(window for _, window, _ in settings)
    if fold_region is None:
        fold_region = (
            max(0, global_start - widest_flank),
            min(len(mrna), global_start + sense_length + widest_flank),
        )
    fold_start, fold_end = fold_region
    fold_sequence = dna_to_rna(mrna[fold_start:fold_end])

    out = {}
    for flank_size, window_size, step in settings:
        cut_start = max(0, global_start - flank_size)
        cut_end = min(len(mrna), global_start + sense_length + flank_size)
        sub_sequence = mrna[cut_start:cut_end]
        if len(sub_sequence) < window_size:
            out[(flank_size, window_size, step)] = np.nan
            continue
        step_to_value = calculate_avg_mfe_per_step(
            sequence=sub_sequence,
            sense_start_in_flank=(global_start - cut_start),
            sense_length=sense_length,
            window_size=window_size,
            steps=[step],
            fold=SharedFold(fold_sequence, cut_start - fold_start, widest_window),
        )
        out[(flank_size, window_size, step)] = step_to_value[step]
    return out
