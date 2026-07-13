import functools
import gc
import logging
import os

import psutil
from numba import njit

logger = logging.getLogger(__name__)


def get_longer_string(s1: str, s2: str) -> str:
    return s1 if len(s1) >= len(s2) else s2


def rna_to_dna(seq) -> str:
    """Uppercase and map the RNA alphabet to DNA (U->T)."""
    return str(seq).upper().replace("U", "T")


def dna_to_rna(seq) -> str:
    """Uppercase and map the DNA alphabet to RNA (T->U)."""
    return str(seq).upper().replace("T", "U")


def _norm_rna_to_dna(seq: str) -> str:
    """Normalize RNA to DNA alphabet (U->T), uppercase, strip whitespace."""
    return rna_to_dna(seq).replace(" ", "").replace("\t", "").replace("\n", "")


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
    return _norm_rna_to_dna(s)


# Legacy dict kept for get_nucleotide_watson_crick (used by numba-compiled code)
WATSON_CRICK_MAP = {"A": "T", "G": "C", "C": "G", "T": "A", "U": "A"}

# Reverse-complement tables (str.translate ~10× faster than per-char dict).
_RC_DNA_TABLE = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")  # output: T
_RC_RNA_TABLE = str.maketrans("ACGTUacgtu", "UGCAAugcaa")  # output: U


def get_antisense(sense: str) -> str:
    """Reverse complement in the DNA alphabet (output uses T)."""
    return sense[::-1].translate(_RC_DNA_TABLE)


def get_antisense_rna(sense: str) -> str:
    """Reverse complement in the RNA alphabet (output uses U, not T)."""
    return sense[::-1].translate(_RC_RNA_TABLE)


def aso_target_rna(aso: str) -> str:
    """The RNA target of an ASO: its reverse complement read as RNA (uppercase)."""
    return get_antisense_rna(aso.upper())


ZERO_CELSIUS_IN_KELVIN = 273.15
# Physiological reference temperature for the nearest-neighbour dG calculations.
BODY_TEMPERATURE_C = 37.0


def celsius_to_kelvin(temp_c: float) -> float:
    return temp_c + ZERO_CELSIUS_IN_KELVIN


def kelvin_to_celsius(temp_k: float) -> float:
    return temp_k - ZERO_CELSIUS_IN_KELVIN


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
