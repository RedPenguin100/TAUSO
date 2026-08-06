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


# The accepted alphabet after U->T normalization. Written down once so the scalar and the
# vectorized validators cannot drift apart.
DNA_BASES = "ACGT"
DNA_BASE_SET = frozenset(DNA_BASES)


def normalize_dna(seq) -> str:
    """Uppercase and map U->T, then check only DNA codes are present -- no N or IUPAC
    ambiguity codes. Returns the normalized sequence. Whitespace is rejected, not stripped."""
    # None/NaN/numbers stringify into plausible-looking letters ("nan" -> "NAN"), which would be
    # reported as a bad base rather than as the wrong type. Sequence-like objects that stringify
    # to the sequence itself (Bio.Seq, numpy.str_) still pass through.
    if seq is None or isinstance(seq, (int, float)):
        raise TypeError(f"sequence must be a string, got {type(seq).__name__}")
    s = str(seq)
    if any(c.isspace() for c in s):
        raise ValueError("sequence must not contain whitespace")
    s = rna_to_dna(s)
    if not s:
        raise ValueError("sequence is empty")
    bad = sorted(set(s) - DNA_BASE_SET)
    if bad:
        raise ValueError(f"sequence has non-ACGU characters: {bad[:5]}")
    return s


def validate_dna_sequences(seq_series) -> None:
    """Raise if any sequence in `seq_series` holds anything but DNA codes -- no N or IUPAC
    ambiguity codes. Expects a Series already uppercased and mapped U->T."""
    invalid = seq_series[seq_series.str.contains(rf"[^{DNA_BASES}]", regex=True, na=True)]
    if len(invalid) > 0:
        raise ValueError(
            f"{len(invalid)} sequence(s) contain non-ACGT characters after upper/U->T normalization "
            f"(only A/C/G/T/U accepted). Examples: {invalid.head(5).tolist()}"
        )


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
