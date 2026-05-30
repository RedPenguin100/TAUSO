"""Sequence-derived liability / toxicity features for antisense oligonucleotides.

Each feature is tied to a distinct, literature-grounded liability mechanism and cited inline.

IMPORTANT: toxicity is not inefficacy. A sequence can be hepatotoxic or immunostimulatory
yet still knock down its target efficiently in a cell line, so these features are not
expected to correlate with Inhibition(%) in general. The exception is the G-quadruplex /
G-run family, which has a direct mechanistic route to *reduced* knockdown (the oligo
aggregates / self-folds instead of hybridising).
"""

import re

from ...common.modifications import get_longest_dna_gap


def _to_dna(seq: str) -> str:
    return seq.upper().replace("U", "T")


def _count_overlapping(s: str, motif: str) -> int:
    return sum(1 for i in range(len(s) - len(motif) + 1) if s[i : i + len(motif)] == motif)


# --- Hepatotoxicity (RNase-H gapmer, target-independent) -----------------------------------
# Burdick AD, Sciabola S, Mantena SR, et al. "Sequence motifs associated with hepatotoxicity
# of locked nucleic acid-modified antisense oligonucleotides." Nucleic Acids Res.
# 2014;42(8):4882-4891. doi:10.1093/nar/gku142. A random-forest model found two trinucleotides,
# TCC and TGC, present only in hepatotoxic sequences. Counted within the deoxy gap here because
# the liability is driven by RNase-H cleavage, which acts on the DNA gap.
_HEPATOTOX_MOTIFS = ("TCC", "TGC")


def hepatotox_gap_motif_count(seq: str, chemical_pattern: str) -> int:
    """Count of hepatotoxicity-associated 3-mers (TCC, TGC) within the longest DNA gap.

    Source: Burdick et al. 2014, Nucleic Acids Res. 42(8):4882, doi:10.1093/nar/gku142.
    """
    if not isinstance(seq, str) or not isinstance(chemical_pattern, str) or len(seq) != len(chemical_pattern):
        return 0
    start, end, length = get_longest_dna_gap(chemical_pattern)
    if length < 3:
        return 0
    gap = _to_dna(seq[start:end])
    return sum(_count_overlapping(gap, m) for m in _HEPATOTOX_MOTIFS)


def hepatotox_motif_count(seq: str) -> int:
    """Count of hepatotoxicity-associated 3-mers (TCC, TGC) over the whole ASO.

    Source: Burdick et al. 2014, Nucleic Acids Res. 42(8):4882, doi:10.1093/nar/gku142
    (the paper scored these as sequence-level motifs).
    """
    s = _to_dna(seq)
    return sum(_count_overlapping(s, m) for m in _HEPATOTOX_MOTIFS)


# --- Immunostimulation via TLR9 ------------------------------------------------------------
# Krieg AM, Yi AK, Matson S, et al. "CpG motifs in bacterial DNA trigger direct B-cell
# activation." Nature 1995;374(6522):546-549. doi:10.1038/374546a0. Identified unmethylated
# CpG in an RR-CG-YY context as the broad stimulatory consensus in mouse B cells; the four
# most stimulatory hexamers were GACGTC, GACGTT, AACGTC, AACGTT, vs the least stimulatory
# TCCGGA, TCCGGG, CCCGGA, CCCGGG.
#
# Hartmann G, Krieg AM. "Mechanism and function of a newly identified CpG DNA motif in
# human primary B cells." J Immunol 2000;164(2):944-953. doi:10.4049/jimmunol.164.2.944.
# Identified GTCGTT (and TTCGTT, AACGTT) as the human-optimal CpG-B hexamer; CpG ODN 2006
# (agatolimod) carries three GTCGTT motifs.
_TLR9_HUMAN_HEXAMER = re.compile(r"GTCGTT")


def cpg_count(seq: str) -> int:
    """Number of CpG dinucleotides (TLR9 immunostimulation proxy).

    Source: Krieg et al. 1995, Nature 374(6522):546, doi:10.1038/374546a0.
    """
    s = _to_dna(seq)
    return sum(1 for i in range(len(s) - 1) if s[i : i + 2] == "CG")


def tlr9_human_motif_count(seq: str) -> int:
    """Count of the human-optimal CpG-B hexamer GTCGTT.

    Source: Hartmann & Krieg 2000, J Immunol 164(2):944, doi:10.4049/jimmunol.164.2.944.
    """
    return len(_TLR9_HUMAN_HEXAMER.findall(_to_dna(seq)))


# --- G-quadruplex / G-aggregation (the only efficacy-relevant liability here) --------------
# Bedrat A, Lacroix L, Mergny JL. "Re-evaluation of G-quadruplex propensity with G4Hunter."
# Nucleic Acids Res. 2016;44(4):1746-1759. doi:10.1093/nar/gkw006. Per-base score: a run of
# k G's scores +min(k,4) per base, a run of k C's scores -min(k,4); the G4Hunter value is the
# mean over a window. Positive = G-rich / quadruplex-prone.
def _g4hunter_scores(seq: str) -> list[int]:
    s = _to_dna(seq)
    n = len(s)
    scores = [0] * n
    i = 0
    while i < n:
        b = s[i]
        if b in "GC":
            j = i
            while j < n and s[j] == b:
                j += 1
            scores[i:j] = [min(j - i, 4) * (1 if b == "G" else -1)] * (j - i)
            i = j
        else:
            i += 1
    return scores


def g4hunter_mean(seq: str) -> float:
    """G4Hunter score over the whole ASO (Bedrat et al. 2016, NAR 44(4):1746, doi:10.1093/nar/gkw006)."""
    scores = _g4hunter_scores(seq)
    return sum(scores) / len(scores) if scores else 0.0


def g4hunter_max(seq: str, window: int = 10) -> float:
    """Maximum windowed G4Hunter score (most G-quadruplex-prone stretch).

    Source: Bedrat et al. 2016, Nucleic Acids Res. 44(4):1746, doi:10.1093/nar/gkw006.
    """
    scores = _g4hunter_scores(seq)
    n = len(scores)
    if n == 0:
        return 0.0
    w = min(window, n)
    return max(sum(scores[k : k + w]) / w for k in range(n - w + 1))


def max_g_run(seq: str) -> int:
    """Length of the longest run of consecutive G's (G-aggregation proxy)."""
    runs = re.findall(r"G+", _to_dna(seq))
    return max((len(r) for r in runs), default=0)


def gggg_motif_count(seq: str) -> int:
    """Overlapping count of GGGG occurrences (literal 4+ G-tract motif).

    Complement to :func:`g4hunter_max`: G4Hunter is a windowed mean, so an
    isolated GGGG that doesn't dominate the window can score low even though
    a 4-G tract is biologically a strong G-quadruplex / G-aggregation flag on
    its own. This is a literal motif count, reported overlapping (GGGG = 1,
    GGGGG = 2, GGGGGG = 3) so longer G-runs score proportionally higher.
    """
    return _count_overlapping(_to_dna(seq), "GGGG")


# --- Acute CNS / neurotoxicity (non-hybridisation) -----------------------------------------
# Hagedorn PH, Brown JM, Easton A, et al. "Acute Neurotoxicity of Antisense Oligonucleotides
# After Intracerebroventricular Injection Into Mouse Brain Can Be Predicted from Sequence
# Features." Nucleic Acid Ther. 2022;32(3):151-162. doi:10.1089/nat.2021.0071. Guanines near
# the 3' end and a short 3'-terminal G-free stretch are associated with acute neurotoxicity.
def three_prime_g_fraction(seq: str, window: int = 4) -> float:
    """Fraction of G in the 3'-terminal `window` bases.

    Source: Hagedorn et al. 2022, Nucleic Acid Ther. 32(3):151, doi:10.1089/nat.2021.0071.
    """
    s = _to_dna(seq)
    if len(s) < window:
        window = len(s)
    if window == 0:
        return 0.0
    tail = s[-window:]
    return tail.count("G") / window


def three_prime_g_free_length(seq: str) -> int:
    """Length of the 3'-terminal stretch containing no G (short stretch = neurotox-associated).

    Source: Hagedorn et al. 2022, Nucleic Acid Ther. 32(3):151, doi:10.1089/nat.2021.0071.
    """
    s = _to_dna(seq)
    n = 0
    for base in reversed(s):
        if base == "G":
            break
        n += 1
    return n
