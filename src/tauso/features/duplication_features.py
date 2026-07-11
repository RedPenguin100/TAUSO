"""Target-duplication features."""

from ..util import get_antisense_rna


def aso_target_rna(aso: str) -> str:
    """The RNA target of an ASO: its reverse complement read as RNA."""
    return get_antisense_rna(aso.upper())


def _exact_starts(mrna: str, target: str, cap: int = 3000) -> list[int]:
    """Start indices of every exact occurrence of ``target`` in ``mrna``."""
    starts = []
    i = mrna.find(target)
    while i != -1 and len(starts) < cap:
        starts.append(i)
        i = mrna.find(target, i + 1)
    return starts


def _approx_starts(mrna: str, target: str, max_mm: int, cap: int = 6000) -> list[int]:
    """Start indices of occurrences of ``target`` in ``mrna`` within ``max_mm`` mismatches.

    Splits ``target`` into ``max_mm + 1`` non-overlapping seeds. By the pigeonhole principle an
    alignment with at most ``max_mm`` mismatches leaves at least one seed exact, so scanning the seeds
    yields every candidate start, which is then verified by counting mismatches.
    """
    L = len(target)
    k = max_mm + 1
    cand = set()
    for c in range(k):
        a, b = c * L // k, (c + 1) * L // k
        seed = target[a:b]
        if len(seed) < 4:
            continue
        i = mrna.find(seed)
        while i != -1 and len(cand) < cap:
            start = i - a
            if 0 <= start <= len(mrna) - L:
                cand.add(start)
            i = mrna.find(seed, i + 1)
    starts = [st for st in cand if sum(x != y for x, y in zip(mrna[st : st + L], target)) <= max_mm]
    starts.sort()
    return starts


def _count_clusters(starts: list[int], min_gap: int) -> int:
    """Number of clusters among sorted ``starts``, merging any two hits less than ``min_gap`` apart."""
    clusters = 0
    last = -(10**9)
    for st in starts:
        if st - last >= min_gap:
            clusters += 1
            last = st
    return clusters


def distinct_matches(mrna: str, target: str, max_mm: int = 1) -> int:
    """Number of distinct occurrences of ``target`` in ``mrna`` within ``max_mm`` mismatches.

    Overlapping or tandem hits closer than one target length apart count once.
    """
    if max_mm == 0:
        starts = _exact_starts(mrna, target)
    else:
        starts = _approx_starts(mrna, target, max_mm)
    return _count_clusters(starts, len(target))
