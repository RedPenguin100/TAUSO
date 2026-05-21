from numba import njit


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
