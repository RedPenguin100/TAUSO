"""Target-duplication features."""

import numpy as np

from ..util import aso_target_rna


def _exact_starts(mrna: str, target: str) -> list[int]:
    """Start indices of every exact occurrence of ``target`` in ``mrna``."""
    starts = []
    i = mrna.find(target)
    while i != -1:
        starts.append(i)
        i = mrna.find(target, i + 1)
    return starts


def _mismatches(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def _one_mismatch_starts(mrna: str, target: str) -> list[int]:
    """Start indices where ``target`` occurs in ``mrna`` with at most one mismatch.

    A single mismatch leaves one half of ``target`` exact, so every hit is anchored by an exact match
    of the left or right half; each anchor is then verified.
    """
    L = len(target)
    half = L // 2
    starts = set()
    for offset, seed in ((0, target[:half]), (half, target[half:])):
        i = mrna.find(seed)
        while i != -1:
            start = i - offset
            if 0 <= start <= len(mrna) - L and _mismatches(mrna[start : start + L], target) <= 1:
                starts.add(start)
            i = mrna.find(seed, i + 1)
    return sorted(starts)


def _count_clusters(starts: list[int], min_gap: int) -> int:
    """Number of clusters among sorted ``starts``, merging any two hits less than ``min_gap`` apart."""
    clusters = 0
    last = None
    for st in starts:
        if last is None or st - last >= min_gap:
            clusters += 1
            last = st
    return clusters


def distinct_matches(mrna: str, target: str, max_mm: int = 1) -> int:
    """Number of distinct occurrences of ``target`` in ``mrna`` within ``max_mm`` (0 or 1) mismatches.

    Overlapping or tandem hits closer than one target length apart count once.
    """
    starts = _exact_starts(mrna, target) if max_mm == 0 else _one_mismatch_starts(mrna, target)
    return _count_clusters(starts, len(target))


def compute_duplications(starts, genes, asos, gene_mrna):
    """Per-ASO counts of near-full-length target copies in the ASO's own pre-mRNA.

    ``gene_mrna`` maps gene name -> uppercase RNA pre-mRNA. ``starts`` are the target positions,
    ``-1`` where the target was not located. Returns (dup_exact, dup_near) arrays (exact and
    <=1-mismatch occurrences), NaN for an unlocated ASO whose gene sequence is also missing. A
    located ASO whose gene has no pre-mRNA is an inconsistency and raises.
    """
    n = len(genes)
    dup_exact = np.full(n, np.nan)
    dup_near = np.full(n, np.nan)
    for i in range(n):
        mrna = gene_mrna.get(genes[i])
        if mrna is None:
            if starts[i] < 0:
                continue
            raise KeyError(f"{genes[i]}: located ASO (sense_start {starts[i]}) but pre-mRNA is missing")
        target = aso_target_rna(asos[i])
        dup_exact[i] = distinct_matches(mrna, target, 0)
        dup_near[i] = distinct_matches(mrna, target, 1)
    return dup_exact, dup_near
