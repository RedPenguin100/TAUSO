"""Base composition of the pre-mRNA flanks around the ASO target site."""

import numpy as np

from ..util import _to_str_seq
from .sequence.seq_features import at_skew, get_gc_content

FLANK = 20


def compute_flank_composition(starts, lengths, genes, gene_to_data, flank=FLANK):
    """AT-skew and GC-content of the +/-``flank`` nt pre-mRNA flanks around each ASO target site.

    ``starts`` are the target positions in each gene's pre-mRNA, ``-1`` where the target was not
    located. Returns (at_skew, gc_content) arrays, NaN for the unlocated ASOs. A located ASO whose
    gene has no sequence is an inconsistency and raises.
    """
    n = len(starts)
    at_skew_out = np.full(n, np.nan)
    gc_out = np.full(n, np.nan)
    seq_by_gene = {}
    for i in range(n):
        start = starts[i]
        if start < 0:
            continue
        gene = genes[i]
        if gene not in seq_by_gene:
            locus = gene_to_data.get(gene)
            if locus is None:
                raise KeyError(f"{gene}: located ASO (sense_start {start}) but gene sequence is missing")
            seq_by_gene[gene] = _to_str_seq(locus.full_mrna)
        pre = seq_by_gene[gene]
        end = start + lengths[i]
        seq = pre[max(0, start - flank) : start] + pre[end : end + flank]
        at_skew_out[i] = at_skew(seq)
        gc_out[i] = get_gc_content(seq)
    return at_skew_out, gc_out
