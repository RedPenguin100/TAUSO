from typing import Iterable
from collections import defaultdict
from itertools import repeat
import math
from multiprocess.pool import Pool
import numpy as np

# Constants
pref_type = np.uint32  # sequence length of up to 4294967295
ref_index_type = np.uint16  # no. of sequences of up to 65535
########################
# Utility function
def is_str_iter(obj):
    return not isinstance(obj, str) and isinstance(obj, Iterable)
########################

 # Suffix array building
def build_suffix_array(ref, pos_spec=True, n_jobs=None):
    """
    Returns a suffix array dict with the following keys:

        ref: list of strings
        ind: index in ref
        pos: position (from start) in string

        Additional fields for position-specific Chimera:
        pos_from_stop: position (from stop) in string
        sorted_start: indices of suffixes sorted by `pos`
        sorted_stop: indices of suffixes sorted by `pos_from_stop`
    """
    if not is_str_iter(ref):
        ref = [ref]

    with Pool(n_jobs) as pool:
        # SA per reference sequence
        SA = pool.starmap(build_single_suffix_array, enumerate(ref))
        # Aggregate, merge-sort style
        while len(SA) > 1:
            merged_SA = pool.starmap(merge_arrays, zip(SA[0::2], SA[1::2], repeat(ref)))
            if len(merged_SA) < len(SA) / 2:
                merged_SA.append(SA[-1])
            SA = merged_SA

    SA = {
        'ref': ref,
        'ind': SA[0][1].astype(ref_index_type),
        'pos': SA[0][0].astype(pref_type),
        'homologs': set()
    }

    if not pos_spec:
        return SA

    # Additional fields for position-specific Chimera
    lens = np.array([len(r) for r in ref])
    SA['pos_from_stop'] = SA['pos'] - lens[SA['ind']]

    return SA
##################

# Functions related to suffix arrays
def longest_prefix(key, SA, max_len=40):
    """ Standard LCS search with a homolog sequence filter based on max_len. """
    max_len = min([len(key), max_len])
    where = search_suffix(key, SA)
    where = min([max([0, where]), SA['pos'].size - 1])
    n = np.inf  # Prefix length

    while n > max_len:
        if np.isfinite(n):
            SA['homologs'].add(SA['ind'][pind])

        nei = get_neighbors(SA, where)  # Adjacent suffixes
        if not len(nei):
            pind = -1
            pref = ''
            return pref, pind

        nei_count = [count_common(key, SA, n) for n in nei]
        nei_max = np.argmax(nei_count)
        pind = nei[nei_max]  # Prefix index in SA
        n = nei_count[nei_max]

    pref = key[:nei_count[nei_max]]  # Prefix string
    return pref, pind
####################
def search_suffix(key, SA, top=None, bottom=None):
    """ Using binary search to find the position where key should be inserted. """
    if top is None:
        top = 0
    if bottom is None:
        bottom = SA['pos'].size - 1
    while top < bottom:
        mid = (top + bottom) // 2
        if is_key_greater_than(key, SA, mid):
            top = mid + 1
        else:
            bottom = mid

    nS = SA['pos'].size - 1
    if top < nS:
        top = skip_masked_suffix(SA, top, +1)
    if top == nS and is_suffix_masked(SA, top):
        top = skip_masked_suffix(SA, top, -1)
    if is_key_greater_than(key, SA, top):
        top += 1

    return top
###################
def get_neighbors(SA, ind):
    """ Returns neighbors of a suffix in SA. """
    nS = SA['pos'].size
    lo = skip_masked_suffix(SA, max([0, ind - 1]), -1)
    hi = skip_masked_suffix(SA, min([nS - 1, ind + 1]), +1)
    neis = [n for n in [lo, ind, hi] if not is_suffix_masked(SA, n)]
    return neis
#####################
def count_common(key, SA, i):
    """ Counts the number of common characters between a key and a suffix. """
    suf = get_suffix(SA, i)[:len(key)]
    key = key[:len(suf)]
    if key == suf:
        return len(key)
    for j, (a, b) in enumerate(zip(key, suf)):
        if a != b:
            break
    return j
###################

def is_key_greater_than(key, SA, i):
    return get_suffix(SA, i) < key
#####################
def get_suffix(SA, i):
    return SA['ref'][SA['ind'][i]][SA['pos'][i]:]
#####################
def skip_masked_suffix(SA, ind, step, min_ind=0, max_ind=None):
    """ Skips masked suffixes in SA. """
    if max_ind is None:
        max_ind = SA['pos'].size

    while is_suffix_masked(SA, ind) and (min_ind <= ind + step) and (ind + step < max_ind):
        ind = ind + step

    return ind
#####################
def is_suffix_masked(SA, ind):
    if 'win_start' not in SA and 'win_stop' not in SA:
        return SA['ind'][ind] in SA['homologs']

    in_win = False
    if 'win_start' in SA:
        in_win |= SA['pos'][ind] in SA['win_start']
    if 'win_stop' in SA:
        in_win |= SA['pos_from_stop'][ind] in SA['win_stop']

    return (not in_win) or (SA['ind'][ind] in SA['homologs'])
#########################
def merge_arrays(SA1, SA2, ref):
    """ Merges two suffix arrays into one. """
    i1 = 0
    i2 = 0
    insert = 0

    n1 = SA1.shape[1]
    n2 = SA2.shape[1]
    if n1 == 0:
        return SA2
    if n2 == 0:
        return SA1
    outSA = np.zeros((2, n1 + n2), dtype=pref_type)

    suf1 = get_raw_suffix(SA1, i1, ref)
    suf2 = get_raw_suffix(SA2, i2, ref)
    while (i1 < n1) and (i2 < n2):
        is_greater = suf1 < suf2
        if is_greater:
            outSA[:, insert] = SA1[:, i1]
            i1 += 1
            if i1 < n1:
                suf1 = get_raw_suffix(SA1, i1, ref)
        else:
            outSA[:, insert] = SA2[:, i2]
            i2 += 1
            if i2 < n2:
                suf2 = get_raw_suffix(SA2, i2, ref)
        insert += 1

    if i1 < n1:
        outSA[:, insert:] = SA1[:, i1:]
    if i2 < n2:
        outSA[:, insert:] = SA2[:, i2:]

    return outSA
#########################
def get_raw_suffix(SA, i, ref):
    return ref[SA[1, i]][SA[0, i]:]
######################
def build_single_suffix_array(i, string):
    return np.vstack([build_suffix_array_ManberMyers(string), len(string) * [i]])
#####################
# Suffix array using Manber-Myers algorithm
def build_suffix_array_ManberMyers(string):
    result = []

    def sort_bucket(string, bucket, order=1):
        d = defaultdict(list)
        for i in bucket:
            key = string[i:i + order]
            d[key].append(i)
        for k, v in sorted(d.items()):
            if len(v) > 1:
                sort_bucket(string, v, order * 2)
            else:
                result.append(v[0])
        return result

    return sort_bucket(string, (i for i in range(len(string))))
####################


def calc_suffix_array(reference_set):
    suffix_array = build_suffix_array(reference_set)
    return suffix_array


