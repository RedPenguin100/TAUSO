import numpy as np
import ViennaRNA as RNA

from ...common.thermo_consts import RT
from ...util import dna_to_rna


def _unpaired_table(sequence, u_max, max_bp_span):
    """
    Return `up[j][u]` - the probability the `u` nucleotides ending at `j` are unpaired.
    """
    span = len(sequence) if max_bp_span is None else max_bp_span
    return RNA.pfl_fold_up(sequence, u_max, len(sequence), span)


def reduce_energies(energies, reducer):
    """Collapse the windows an anchoring selected to one number.

    `mean` is how open the target is on average. `std` is how uneven that is across it:
    two sites can share a mean while one is uniformly half-open and the other alternates
    between a closed stretch and an open one, and only the spread separates them. A
    single-window anchoring has no spread, so `std` over one window is always zero.
    """
    if reducer == "mean":
        return float(np.mean(energies))
    if reducer == "std":
        return float(np.std(energies))
    raise ValueError(f"unknown reducer {reducer!r}; expected one of mean, std")


def window_starts(anchor, sense_start, sense_length, open_len):
    """0-based starts of the windows an anchoring averages over.

    `a5` and `a3` sweep the whole target, one window per position: `a5` takes those
    starting inside it, so they hang over its 3' edge; `a3` takes those ending inside it,
    hanging over the 5' edge instead. A window wider than the target cannot sit inside it,
    so which edge it overhangs is a choice, and the two see different flanking context.

    `aso5end` and `aso3end` take a single window flush with one end of the target. The
    oligo binds antiparallel, so its 5' terminus pairs with the target's 3' side and its
    3' terminus with the target's 5' side; RNase H1 then cleaves in a fixed register from
    the ASO-3' end, which is why the two ends are not interchangeable.
    """
    sense_end = sense_start + sense_length
    if anchor == "a5":
        return range(sense_start, sense_end)
    if anchor == "a3":
        return range(sense_start - open_len + 1, sense_end - open_len + 1)
    if anchor == "aso3end":
        return range(sense_start, sense_start + 1)
    if anchor == "aso5end":
        return range(sense_end - open_len, sense_end - open_len + 1)
    raise ValueError(f"unknown anchor {anchor!r}; expected one of a5, a3, aso5end, aso3end")


def opening_energies(unpaired, starts, open_len):
    """Opening energy at each window start.

    A window starting at 0-based `i` ends at 1-based `i + open_len`. Starts running off
    either end of the cut have no entry, as do those whose probability underflows.
    """
    out = []
    for start in starts:
        end = start + open_len
        if start < 0 or end >= len(unpaired):
            continue
        probability = unpaired[end][open_len]
        if probability and probability > 0.0:
            out.append(-RT * np.log(probability))
    return out


def calculate_avg_access_per_setting(mrna, global_start, sense_length, settings):
    """Accessibility of the target, one value per (flank, max_bp_span, open_len, anchor, reducer).

    Accessibility is the free energy needed to hold a stretch of the target
    single-stranded. Unlike MFE it is an ensemble quantity: it sums over every
    structure the cut can adopt, weighted by Boltzmann probability.

    Settings sharing a (flank, max_bp_span) fold the same cut, so they are grouped and
    folded once at the widest opening length in the group; every opening length and
    anchoring is then a free readout of that one fold. The value is the mean opening
    energy over the anchoring's windows.

    Returns a dict keyed by the setting tuples, NaN where the cut is too short to hold
    a window or every window falls outside it.
    """
    by_cut = {}
    for flank, max_bp_span, open_len, anchor, reducer in settings:
        by_cut.setdefault((flank, max_bp_span), []).append((open_len, anchor, reducer))

    out = {}
    for (flank, max_bp_span), reads in by_cut.items():
        cut_start = max(0, global_start - flank)
        cut_end = min(len(mrna), global_start + sense_length + flank)
        cut = dna_to_rna(mrna[cut_start:cut_end])
        u_max = max(open_len for open_len, _, _ in reads)
        if len(cut) < u_max:
            for open_len, anchor, reducer in reads:
                out[(flank, max_bp_span, open_len, anchor, reducer)] = np.nan
            continue

        unpaired = _unpaired_table(cut, u_max, max_bp_span)
        sense_start = global_start - cut_start
        for open_len, anchor, reducer in reads:
            starts = window_starts(anchor, sense_start, sense_length, open_len)
            energies = opening_energies(unpaired, starts, open_len)
            value = reduce_energies(energies, reducer) if energies else np.nan
            out[(flank, max_bp_span, open_len, anchor, reducer)] = value
    return out
