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
        if not hasattr(cls, "_str_map"):
            cls._str_map = {
                "unknown": cls.UNKNOWN,
                "protein_coding": cls.PROTEIN_CODING,
                "lncrna": cls.LNCRNA,  # Lowercase for safe matching
                "mirna": cls.MIRNA,
                "snorna": cls.SNO_RNA,
                "non_coding": cls.NON_CODING,  # FIXED: Was mapped to SNO_RNA
                "unannotated": cls.UNANNOTATED,
            }
        # Normalize to lowercase so exact casing from the GFF doesn't break the mapping
        clean_str = type_str.lower() if type_str else ""
        return cls._str_map.get(clean_str, cls.UNKNOWN)


class StrandType(IntEnum):
    POS = 1
    NEG = 2

    @classmethod
    def from_string(cls, type_str):
        _map = {"+": cls.POS, "-": cls.NEG}
        return _map[type_str]


class LocusInfo:
    __slots__ = [
        "exon_indices",
        "intron_indices",
        "stop_codons",
        "five_prime_utr",
        "three_prime_utr",
        "full_mrna",
        "gene_start",
        "gene_end",
        "strand",
        "gene_type",
        "utr_indices",
    ]

    def __init__(self, seq=None):
        self.exon_indices = []
        self.intron_indices = []
        self.stop_codons = []
        self.five_prime_utr = ""
        self.three_prime_utr = ""
        self.full_mrna = None
        self.gene_start = None
        self.gene_end = None
        self.strand = None
        self.gene_type = None
        self.utr_indices = []

        if seq is not None:
            self.exon_indices = [(0, len(seq))]
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

    @property
    def exons(self):
        return [self._get_sequence_slice(start, end) for start, end in self.exon_indices]

    @property
    def introns(self):
        return [self._get_sequence_slice(start, end) for start, end in self.intron_indices]

    @property
    def exon_concat(self):
        return "".join(self.exons)

    def __repr__(self):
        # Since __dict__ no longer exists with __slots__,
        # we manually iterate the slots for the representation.
        lines = ["LocusInfo:"]
        for field in self.__slots__:
            lines.append(f"  {field}: {getattr(self, field)}")
        return "\n".join(lines)