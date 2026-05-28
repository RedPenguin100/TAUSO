from .seq_features import get_gc_content

# High-affinity sugar wings of a gapmer: 'M' = 2'-MOE, 'C' = cEt, 'L' = LNA.
_WING_MARKERS = "MCL"


def gap_gc_content(seq: str, pattern: str) -> float:
    """Calculates GC content strictly within the DNA gap ('d')."""
    gap_bases = "".join(base for base, pat in zip(seq, pattern) if pat == "d")
    return get_gc_content(gap_bases)


def wing_gap_gc_delta(seq: str, pattern: str) -> float:
    """
    Calculates the difference in GC content between the high-affinity wings
    ('M' = 2'-MOE, 'C' = cEt, 'L' = LNA) and the central DNA gap ('d').
    """
    wing_bases = "".join(base for base, pat in zip(seq, pattern) if pat in _WING_MARKERS)
    gap_bases = "".join(base for base, pat in zip(seq, pattern) if pat == "d")

    if not wing_bases or not gap_bases:
        return 0.0

    return get_gc_content(wing_bases) - get_gc_content(gap_bases)
