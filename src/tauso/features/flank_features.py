"""Base composition of the pre-mRNA flanks around the ASO target site."""

import numpy as np

FLANK = 20


def _composition(flank: str) -> tuple[float, float]:
    """AT-skew (A-T)/(A+T) and GC-content (G+C)/(A+C+G+T) of a flank sequence."""
    a, t, g, c = flank.count("A"), flank.count("T"), flank.count("G"), flank.count("C")
    at_skew = (a - t) / (a + t) if (a + t) else 0.0
    gc = (g + c) / (a + c + g + t) if (a + c + g + t) else 0.0
    return at_skew, gc


def compute_flank_composition(starts, lengths, genes, gene_registry, flank=FLANK):
    """AT-skew and GC-content of the +/-``flank`` nt pre-mRNA flanks around each ASO target site.

    ``gene_registry`` maps gene name -> a record with a ``pre_mrna_sequence`` (the transcript-sense
    pre-mRNA that ``starts`` indexes into). Returns (at_skew, gc_content) arrays, NaN where the
    gene/start is unavailable.
    """
    n = len(starts)
    at_skew = np.full(n, np.nan)
    gc = np.full(n, np.nan)
    for i in range(n):
        idx = starts[i]
        if not np.isfinite(idx) or idx < 0:
            continue
        rec = gene_registry.get(genes[i]) if hasattr(gene_registry, "get") else None
        pre = rec.get("pre_mrna_sequence") if isinstance(rec, dict) else None
        if not pre:
            continue
        idx = int(idx)
        ln = int(lengths[i])
        seq = (pre[max(0, idx - flank) : idx] + pre[idx + ln : idx + ln + flank]).upper()
        if seq:
            at_skew[i], gc[i] = _composition(seq)
    return at_skew, gc
