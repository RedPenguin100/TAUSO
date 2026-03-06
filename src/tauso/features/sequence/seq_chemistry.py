def gap_gc_content(seq: str, pattern: str) -> float:
    """Calculates GC content strictly within the DNA gap ('d')."""
    gap_bases = [base.upper() for base, pat in zip(seq, pattern) if pat == 'd']
    if not gap_bases:
        return 0.0
    gc_count = sum(1 for b in gap_bases if b in 'GC')
    return gc_count / len(gap_bases)


def wing_gap_gc_delta(seq: str, pattern: str) -> float:
    """
    Calculates the difference in GC content between the wings ('M', 'C') 
    and the central gap ('d').
    """
    wing_bases = [base.upper() for base, pat in zip(seq, pattern) if pat in 'MC']
    gap_bases = [base.upper() for base, pat in zip(seq, pattern) if pat == 'd']

    if not wing_bases or not gap_bases:
        return 0.0

    wing_gc = sum(1 for b in wing_bases if b in 'GC') / len(wing_bases)
    gap_gc = sum(1 for b in gap_bases if b in 'GC') / len(gap_bases)

    return wing_gc - gap_gc
