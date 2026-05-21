import functools
import gc
import logging
import os

import psutil
from numba import njit

logger = logging.getLogger(__name__)


def get_longer_string(s1: str, s2: str) -> str:
    return s1 if len(s1) >= len(s2) else s2


def _norm_rna_to_dna(seq: str) -> str:
    """Normalize RNA to DNA alphabet (U->T), uppercase, strip whitespace."""
    return str(seq).upper().replace("U", "T").replace(" ", "").replace("\t", "").replace("\n", "")


def _to_str_seq(x) -> str:
    """
    Coerce sequence-like (list/np.array/Series) or string to a clean uppercase DNA string.
    Converts U->T and strips whitespace. Ensures slicing returns a plain string (avoids pandas iterable assignment).
    """
    if isinstance(x, str):
        s = x
    else:
        try:
            s = "".join(list(x))
        except Exception:
            s = str(x)
    return s.replace(" ", "").replace("\t", "").replace("\n", "").replace("U", "T").upper()


# Legacy dict kept for get_nucleotide_watson_crick (used by numba-compiled code)
WATSON_CRICK_MAP = {"A": "T", "G": "C", "C": "G", "T": "A", "U": "A"}

# str.translate tables — ~10× faster than per-char dict lookup
_DNA_TABLE = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")  # output: T
_RC_RNA_TABLE = str.maketrans("ACGTUacgtu", "UGCAAugcaa")  # output: U


def get_antisense(sense: str) -> str:
    """Reverse complement in DNA alphabet (output uses T)."""
    return sense[::-1].translate(_DNA_TABLE)


def get_antisense_rna(sense: str) -> str:
    """Reverse complement in RNA alphabet (output uses U, not T)."""
    return sense[::-1].translate(_RC_RNA_TABLE)


COMP_U_TABLE = bytes.maketrans(b"ACGTUacgtu", b"UGCAaugcaa")


def get_antisense_u(seq: str) -> str:
    return seq.translate(COMP_U_TABLE)[::-1]


@njit
def get_nucleotide_watson_crick(nucleotide):
    if nucleotide == "A":
        return "T"
    if nucleotide == "G":
        return "C"
    if nucleotide == "C":
        return "G"
    if nucleotide == "U" or nucleotide == "T":
        return "A"
    raise ValueError(f"Unknown nucleotide {nucleotide}")


def log_memory_usage(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        process = psutil.Process(os.getpid())
        mem_before = process.memory_info().rss / (1024**2)

        # Force GC to get a clean baseline
        gc.collect()

        result = func(*args, **kwargs)

        gc.collect()
        mem_after = process.memory_info().rss / (1024**2)
        logger.debug(
            f"--- [MEM] {func.__name__} | Before: {mem_before:.1f}MB | After: {mem_after:.1f}MB | Change: {mem_after - mem_before:.1f}MB ---"
        )
        return result

    return wrapper
