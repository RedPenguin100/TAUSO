import math

import numpy as np
import pandas as pd

from Bio.SeqUtils import gc_fraction

from .rna_access import RNAAccess
from ..._raccess.core import find_raccess


def get_sense_with_flanks(pre_mrna: str, sense_start: int, sense_length: int, flank_size: int) -> str:
    """
    Re  turns the sense sequence with `flank_size` nucleotides on each side (if available).
    If near the edge, it will not go out of bounds.

    Parameters:
    - pre_mrna: The full pre-mRNA sequence (5' -> 3')
    - sense_start: Start index of the sense sequence within pre_mrna
    - sense_length: Length of the sense sequence (usually same as antisense length)
    - flank_size: Number of nucleotides to include on each side (upstream and downstream)

    Returns:
    - str: The flanked sense sequence
    """
    # Ensure indices are within bounds
    start = max(0, sense_start - flank_size)
    end = min(len(pre_mrna), sense_start + sense_length + flank_size)

    return pre_mrna[start:end]

def get_cache(seed_sizes, access_size):

    # Mimick original logic
    seed_sizes = list(filter(lambda x: x <= access_size, seed_sizes))


    cache = {}
    for super_seed_size in seed_sizes:
        step = super_seed_size // 2
        rel_offsets = list(range(0, access_size - super_seed_size, step))
        rel_offsets.append(access_size - super_seed_size)

        n_values = math.ceil((access_size - super_seed_size) / step) + 1
        assert len(rel_offsets) == n_values

        weights = np.ones(n_values)
        if n_values > 1:
            last_weight = (rel_offsets[-1] - rel_offsets[-2]) / step
            weights[-1] = last_weight

        cache[super_seed_size] = (np.array(rel_offsets), weights, access_size / super_seed_size)
    return cache




# noinspection DuplicatedCode
class AccessCalculator(object):

    @classmethod
    def calc_gc_info(cls, rna_seq, segment_size):

        trigger_mrna_size = len(rna_seq)

        indexes = list(range(0, trigger_mrna_size - segment_size + 1))
        trigger_segments = [rna_seq[i:i + segment_size] for i in indexes]

        gc_segments = list(map(gc_fraction, trigger_segments))

        d = {
            'trigger_seq': trigger_segments,
            'gc': gc_segments,
        }

        df = pd.DataFrame(d, index=indexes)

        return df

    @classmethod
    def calc_access_energies(
            cls, rna_seq, access_size, seed_sizes, max_span, uuid_str, cache=None):

        rna_size = len(rna_seq)

        seed_sizes = list(filter(lambda x: x <= access_size, seed_sizes))
        assert seed_sizes
        cls.rna_access = RNAAccess(seed_sizes, max_span, find_raccess())
        ra = cls.rna_access
        ra.set_uuid_for_web(uuid_str)
        access_query = [('rna', rna_seq)]
        res = ra.calculate(access_query)

        access_res = res['rna']
        cols_np = {s: access_res[s].to_numpy() for s in seed_sizes}

        ind_info_list = []
        for pos in range(0, rna_size - access_size + 1):
            pos_info = {}
            for super_seed_size in seed_sizes:
                step = super_seed_size // 2
                # n_samples = trigger_binding_site_size / step
                rel_offsets = np.arange(0, access_size - super_seed_size, step)
                rel_offsets = np.append(rel_offsets, access_size - super_seed_size)

                abs_offsets = rel_offsets + pos
                vals = cols_np[super_seed_size][len(cols_np[super_seed_size]) - 1 - abs_offsets]
                bind_energies = pd.Series(vals, index=abs_offsets)

                norm_factor = access_size / super_seed_size
                norm_bind_energies = bind_energies * norm_factor

                # fix last weight relatively to the overlap with the one before
                n_values = math.ceil((access_size - super_seed_size) / step) + 1
                assert len(rel_offsets) == n_values
                weights = [1.0] * n_values
                if len(weights) > 1:
                    last_weight = (rel_offsets[-1] - rel_offsets[-2]) / step
                    weights[-1] = last_weight

                fixed_weight_energies = (np.array(norm_bind_energies) * np.array(weights))
                avg_energy = np.sum(fixed_weight_energies) / np.sum(weights)
                min_energy = np.min(fixed_weight_energies)
                avg_id = f"{super_seed_size}_avg"
                min_id = f"{super_seed_size}_min"
                pos_info.update({avg_id: avg_energy, min_id: min_energy})

            ind_info_list.append((pos, pos_info))

        indexes = list(zip(*ind_info_list))[0]
        records = list(zip(*ind_info_list))[1]
        df = pd.DataFrame(records, index=indexes)
        return df

    @classmethod
    def calc(
            cls, rna_seq, access_size,
            min_gc, max_gc, gc_ranges,
            access_win_size, access_seed_sizes,
            uuid_str=None, temperature=None, cache=None):
        """
        :param rna_seq: target mrna sequence
        :param access_size: target  access size
        :param min_gc: min gc ratio integer between 0 and 100
        :param max_gc: max gc ratio integer between 0 and 100
        :param gc_ranges: if not 1 it should be integer of gc sub ranges between min_gc to max_gc
        :param access_win_size: the sliding window r_access will use for seeking folding interactions
        :param access_seed_sizes: the seed sizes we simulate to check for accessibility segments
        :param uuid_str: RNAAccess module create temporal file so need uuid prefix for parallel run
        :param temperature: temperature for bind energies calculations
        """
        assert len(rna_seq) > 1
        assert len(rna_seq) >= access_size

        rna_seq = rna_seq.upper().replace('T', 'U')

        # gc filter
        gc_info = cls.calc_gc_info(rna_seq, access_size)

        access_energies = cls.calc_access_energies(
            rna_seq, access_size, access_seed_sizes, access_win_size, uuid_str, cache=cache)

        # ae_col1 = f"{access_seed_size}_avg"
        # ae_col2 = f"{access_seed_size * 2}_avg"
        selected_cols = [f"{access_seed_size}_avg" for access_seed_size in access_seed_sizes]
        filtered_access_energies = access_energies.loc[:, access_energies.columns.isin(selected_cols)]

        access_energies['avg_access'] = filtered_access_energies.mean(axis=1)

        assert gc_ranges >= 1
        gc_values = np.linspace(min_gc, max_gc, num=gc_ranges + 1)
        gc_values /= 100

        df_list = []
        for i in range(gc_ranges):
            inc = 'left' if (i + 1 < gc_ranges) else 'both'

            min_range_gc = gc_values[i]
            max_range_gc = gc_values[i + 1]
            gc_indexes = gc_info[gc_info['gc'].between(min_range_gc, max_range_gc, inclusive=inc)].index

            gc_access_energies = access_energies.loc[gc_indexes]

            gc_range = f"gc_range_{min_range_gc}_{max_range_gc}"

            gc_access_energies['gc_range'] = gc_range

            df_list.append(gc_access_energies)

        df = pd.concat(df_list, join='outer', axis=0).fillna(float('nan'))

        return df
