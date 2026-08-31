from dataclasses import dataclass

import numpy as np

from ...util import dna_to_rna
from ._vienna_internal import window_mfe


@dataclass(frozen=True, slots=True)
class SharedFold:
    """
    Class to store the large folded sequence info and fetch subsequence folds.

    `sequence` is the large sequence
    `offset` and `subseq_length` identify the subsequence we which to fold

    `max_bp_span` is an optimization parameter, and should be set to the longest window that will be asked for.
                  Smaller values with cause an error, and bigger values would slow compute.
    """

    sequence: str
    offset: int
    subseq_length: int
    max_bp_span: int

    def __post_init__(self):
        if self.offset < 0 or self.offset + self.subseq_length > len(self.sequence):
            raise ValueError(
                f"sub-sequence [{self.offset}:{self.offset + self.subseq_length}] does not fit "
                f"inside a {len(self.sequence)}nt fold"
            )

    def mfe(self, start: int, window_size: int):
        """
        Equivalent to ViennaRNA.fold(self.sequence[self.offset + start, self.offset + start + window_size])
        """
        return window_mfe(self.sequence, self.offset + start, window_size, self.max_bp_span)

    def per_position_avg_energies(self, window_size, start, stop, step):
        """One sliding-window sweep over the fold's sub-sequence; returns the per-position average MFE.

        Read each window in range(start, stop, step) out of `fold`, spread its mfe/window_size
        evenly over the window's positions, then average per position. Positions no window
        covers are NaN.

        Example:

        subseq_length=12, window_size=5, step=3, range(0, 8, 3)

        pos      0  1  2  3  4  5  6  7  8  9 10 11
        i=0     [--------------]                        mfe/5 added to each of 0..4
        i=3                 [--------------]            mfe/5 added to each of 3..7
        i=6                          [--------------]   mfe/5 added to each of 6..10
                                             ^
        counts   1  1  1  2  2  1  2  2  1  1  1  0  -> NaN, no window reached it

        """
        energy_values = np.zeros(self.subseq_length)
        counts = np.zeros(self.subseq_length)
        for i in range(start, stop, step):
            mfe_per_nt = self.mfe(i, window_size) / window_size
            energy_values[i : i + window_size] += mfe_per_nt
            counts[i : i + window_size] += 1
        return np.divide(
            energy_values,
            counts,
            out=np.full_like(energy_values, np.nan),
            where=counts != 0,
        )

    def sense_profile(self, sense_start_in_flank: int, sense_length: int, window_size: int, step: int):
        """
        Per-position average MFE across the sense region energies.

        pos      0  1  2  3  4  5  6  7  8  9 10 11
        sense                 [==========]                sense_start_in_flank .. sense_end_in_flank -1

        i=2         [-------]                             miss: right edge stops short
        i=3            [-------]                          FIRST  i = sense_start - window + 1
        i=8                           [-------]           LAST   i = sense_end - 1
        i=9                              [-------]        miss: left edge starts past the end
        """
        sense_end_in_flank = sense_start_in_flank + sense_length

        # for when the 5' flank is truncated
        min_i = max(0, sense_start_in_flank - window_size + 1)
        # for when the 3' flank is truncated
        max_i = min(self.subseq_length - window_size, sense_end_in_flank - 1)

        # First multiple of `step` >= min_i. Same grid as starting at 0, just filtered.
        start_i = ((min_i + step - 1) // step) * step
        avg_energies = self.per_position_avg_energies(window_size, start_i, max_i + 1, step)
        return avg_energies[sense_start_in_flank:sense_end_in_flank]


def calculate_avg_mfe(sense_start_in_flank, sense_length, window_size, step, fold):
    """Sliding-window MFE averaged over the sense region. NaN if the target is not inside the cut."""
    if not (
        0 <= sense_start_in_flank < fold.subseq_length and sense_start_in_flank + sense_length <= fold.subseq_length
    ):
        return np.nan
    return np.nanmean(fold.sense_profile(sense_start_in_flank, sense_length, window_size, step))


def calculate_end_mfe(
    full_mrna: str,
    sense_start: int,
    sense_length: int,
    flank_size: int,
    window_size: int,
    step: int,
    end_len: int,
    fold_region,
):
    """Sliding-window MFE over the terminal `end_len` target positions at each ASO end,
    plus `std`: how uneven the target's structure is rather than how much of it there is.

    `fold_region` is the slice of `full_mrna` to fold; it must contain this setting's cut.
    """
    cut_start = max(0, sense_start - flank_size)
    cut_end = min(len(full_mrna), sense_start + sense_length + flank_size)
    cut_len = cut_end - cut_start
    if cut_len < window_size or sense_length < end_len:
        return {"aso5end": np.nan, "aso3end": np.nan, "std": np.nan}
    fold_start, fold_end = fold_region
    fold = SharedFold(dna_to_rna(full_mrna[fold_start:fold_end]), cut_start - fold_start, cut_len, window_size)
    profile = fold.sense_profile(sense_start - cut_start, sense_length, window_size, step)
    # The oligo binds antiparallel, so its 5' end sits over the target's 3' end.
    return {
        "aso5end": float(np.nanmean(profile[sense_length - end_len :])),
        "aso3end": float(np.nanmean(profile[:end_len])),
        "std": float(np.nanstd(profile)),
    }


def calculate_avg_mfe_per_setting(full_mrna, global_start, sense_length, settings, fold_region):
    """Run `calculate_avg_mfe` for every (flank, window, step) in `settings`.

    Each setting sweeps its own flank-padded cut around the target site. Those cuts are all
    substrings of the widest one, so one region covering them is folded once and every
    setting reads its windows out of it. Returns a dict keyed by the setting triples.

    `fold_region` is the (start, end) slice of `full_mrna` to fold, usually one region
    shared with neighbouring target sites. It must contain every setting's cut.
    """
    widest_window = max(window for _, window, _ in settings)
    fold_start, fold_end = fold_region
    fold_sequence = dna_to_rna(full_mrna[fold_start:fold_end])

    out = {}
    for flank_size, window_size, step in settings:
        cut_start = max(0, global_start - flank_size)
        cut_end = min(len(full_mrna), global_start + sense_length + flank_size)
        # The sweep is bounded by this setting's own cut. Bounding it by the folded region
        # instead would let windows run past the cut into a neighbouring site's sequence.
        cut_len = cut_end - cut_start
        if cut_len < window_size:
            out[(flank_size, window_size, step)] = np.nan
            continue
        out[(flank_size, window_size, step)] = calculate_avg_mfe(
            sense_start_in_flank=(global_start - cut_start),
            sense_length=sense_length,
            window_size=window_size,
            step=step,
            # Choosing the widest_window for the max_bp_span is the optimization
            fold=SharedFold(fold_sequence, cut_start - fold_start, cut_len, max_bp_span=widest_window),
        )
    return out
