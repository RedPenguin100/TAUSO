import logging
from enum import IntEnum

logger = logging.getLogger(__name__)


class GeneType(IntEnum):
    UNKNOWN = 0
    PROTEIN_CODING = 1
    LNCRNA = 2
    MIRNA = 3
    SNO_RNA = 4
    NON_CODING = 5
    UNANNOTATED = 6

    @classmethod
    def from_string(cls, type_str):
        if not type_str:
            return cls.UNKNOWN

        # We removed the hasattr() check.
        # Direct lookup first (fast path for exactly matching strings).
        # If it fails, do the .lower() fallback.
        mapped = cls._str_map.get(type_str)
        if mapped is not None:
            return mapped

        return cls._str_map.get(type_str.lower(), cls.UNKNOWN)


# Define the map outside the method so it is only created ONCE at import time.
GeneType._str_map = {
    "unknown": GeneType.UNKNOWN,
    "protein_coding": GeneType.PROTEIN_CODING,
    "lncrna": GeneType.LNCRNA,
    "mirna": GeneType.MIRNA,
    "snorna": GeneType.SNO_RNA,
    "non_coding": GeneType.NON_CODING,
    "unannotated": GeneType.UNANNOTATED,
}


class StrandType(IntEnum):
    POS = 1
    NEG = 2

    @classmethod
    def from_string(cls, type_str):
        # We use a fast dict lookup without recreating the dict
        # Using .get() prevents KeyErrors if corrupted data is passed
        return cls._str_map.get(type_str, cls.POS)  # Defaults to POS, or handle as needed


# Created exactly once.
StrandType._str_map = {"+": StrandType.POS, "-": StrandType.NEG}


class LocusInfo:
    __slots__ = (  # Tuples are slightly faster/safer for __slots__ definition than lists
        "_exon_indices",
        "_intron_indices",
        "stop_codons",
        "five_prime_utr",
        "three_prime_utr",
        "full_mrna",
        "gene_start",
        "gene_end",
        "strand",
        "gene_type",
        "utr_indices",
    )

    def __init__(self, seq=None):
        self._exon_indices = []
        self._intron_indices = []
        self.stop_codons = []
        self.five_prime_utr = ""
        self.three_prime_utr = ""

        self.full_mrna = None
        self.gene_start = None
        self.gene_end = None
        self.strand = None
        self.gene_type = None

        if seq is not None:
            self._exon_indices = [(0, len(seq))]
            self.gene_start = 0
            self.gene_end = len(seq)
            self.strand = StrandType.POS
            self.gene_type = GeneType.UNKNOWN
            self.full_mrna = seq

    def _get_sequence_slice(self, start, end):
        if not self.full_mrna or self.gene_start is None or self.gene_end is None:
            return ""

        if self.strand == StrandType.POS:
            # Relative to the 5' end of the + strand (gene_start)
            rel_start = start - self.gene_start
            rel_end = end - self.gene_start
        else:
            rel_start = self.gene_end - end
            rel_end = self.gene_end - start

        return self.full_mrna[rel_start:rel_end]

    def add_exon_indices(self, start, end):
        self._exon_indices.append((start, end))

    def add_intron_indices(self, start, end):
        self._exon_indices.append((start, end))

    @property
    def exons(self):
        return [self._get_sequence_slice(s, e) for s, e in self._exon_indices]

    @property
    def introns(self):
        return [self._get_sequence_slice(s, e) for s, e in self._intron_indices]

    def __repr__(self):
        lines = [f"  {field}: {getattr(self, field)}" for field in self.__slots__]
        return "LocusInfo:\n" + "\n".join(lines)
