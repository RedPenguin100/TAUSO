"""Target-duplication features.

How many distinct near-full-length copies of an ASO's target sequence appear in its own pre-mRNA.
More copies of the target site mean more substrate for RNase H1 cleavage, which is associated with
stronger knockdown (e.g. repeat-expansion targets such as the C9orf72 GGGGCC repeat).
"""

_COMP = {"A": "U", "U": "A", "G": "C", "C": "G", "N": "N"}


def aso_target_rna(aso: str) -> str:
    """The RNA target of an ASO: its reverse complement read as RNA (T mapped to U)."""
    return "".join(_COMP.get(b, "N") for b in reversed(aso.upper().replace("T", "U")))


def distinct_matches(mrna: str, target: str, max_mm: int = 1) -> int:
    """Number of distinct occurrences of ``target`` in ``mrna`` within ``max_mm`` mismatches.

    Sites closer than one target length apart are collapsed into a single cluster, so overlapping /
    tandem hits count once. The exact case (``max_mm == 0``) scans directly; the mismatch case uses a
    seed-and-extend over ``max_mm + 1`` non-overlapping seeds (pigeonhole: an alignment with at most
    ``max_mm`` mismatches leaves at least one seed exact), then verifies each candidate. Internal caps
    (3000 exact hits, 6000 seed candidates) bound pathological repeat regions.
    """
    L = len(target)
    starts = []
    if max_mm == 0:
        i = mrna.find(target)
        while i != -1 and len(starts) < 3000:
            starts.append(i)
            i = mrna.find(target, i + 1)
    else:
        k = max_mm + 1
        cand = set()
        for c in range(k):
            a, b = c * L // k, (c + 1) * L // k
            seed = target[a:b]
            if len(seed) < 4:
                continue
            i = mrna.find(seed)
            while i != -1 and len(cand) < 6000:
                st = i - a
                if 0 <= st <= len(mrna) - L:
                    cand.add(st)
                i = mrna.find(seed, i + 1)
        for st in cand:
            if sum(x != y for x, y in zip(mrna[st:st + L], target)) <= max_mm:
                starts.append(st)
        starts.sort()
    clusters = 0
    last = -(10**9)
    for st in starts:
        if st - last >= L:
            clusters += 1
            last = st
    return clusters
