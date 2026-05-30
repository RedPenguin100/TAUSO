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

        mapped = cls._str_map.get(type_str)
        if mapped is not None:
            return mapped

        return cls._str_map.get(type_str.lower(), cls.UNKNOWN)


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
    NEG = -1

    @classmethod
    def from_string(cls, type_str):
        return cls._str_map.get(type_str, cls.POS)


StrandType._str_map = {"+": StrandType.POS, "-": StrandType.NEG}


class LocusInfo:
    __slots__ = (
        "_exon_indices",
        "_intron_indices",
        "_5utr_indices",
        "_3utr_indices",
        "stop_codons",
        "all_stop_codons",
        "start_codons",
        "all_start_codons",
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
        self._5utr_indices = []
        self._3utr_indices = []
        self.stop_codons = []
        self.all_stop_codons = []
        self.start_codons = []
        self.all_start_codons = []
        self.five_prime_utr = ""
        self.three_prime_utr = ""
        self.full_mrna = None
        self.gene_start = None
        self.gene_end = None
        self.strand = None
        self.gene_type = None
        self.utr_indices = []

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
            rel_start = start - self.gene_start
            rel_end = end - self.gene_start
        else:
            rel_start = self.gene_end - end
            rel_end = self.gene_end - start
        return self.full_mrna[rel_start:rel_end]

    def add_exon_indices(self, start, end):
        self._exon_indices.append((start, end))

    def add_intron_indices(self, start, end):
        self._intron_indices.append((start, end))

    def add_5utr_indices(self, start, end):
        self._5utr_indices.append((start, end))

    def add_3utr_indices(self, start, end):
        self._3utr_indices.append((start, end))

    @property
    def exons(self):
        return [self._get_sequence_slice(s, e) for s, e in self._exon_indices]

    @property
    def introns(self):
        return [self._get_sequence_slice(s, e) for s, e in self._intron_indices]

    @property
    def sense_utr5(self):
        """
        5' UTR sequences in biological 5'->3' order.
        Indices are already stored in biological order (same sort/reverse
        logic applied in locus_loader), so we just slice.
        """
        return [self._get_sequence_slice(s, e) for s, e in self._5utr_indices]

    @property
    def sense_utr3(self):
        """
        3' UTR sequences in biological 5'->3' order.
        """
        return [self._get_sequence_slice(s, e) for s, e in self._3utr_indices]

    def __repr__(self):
        mrna_len = len(self.full_mrna) if self.full_mrna else 0
        strand_sym = "+" if self.strand == StrandType.POS else "-" if self.strand == StrandType.NEG else "?"
        gene_type_name = self.gene_type.name if self.gene_type is not None else None
        return (
            f"LocusInfo(strand={strand_sym}, type={gene_type_name}, "
            f"coords={self.gene_start}-{self.gene_end}, mrna_len={mrna_len}, "
            f"exons={len(self._exon_indices)}, introns={len(self._intron_indices)})"
        )


class LazyLocusInfo:
    """Drop-in replacement for LocusInfo that does not store the full genomic
    sequence up front.

    Instead it keeps the gene's chromosome + coordinates and fetches `full_mrna`
    from the genome FASTA on first access, memoizing it on the instance. The
    public interface is identical to LocusInfo, so consumers are unchanged; the
    win is that memory scales with the number of genes whose sequence is actually
    read, not the whole genome (the full dict holds ~20k coordinate-only objects).
    """

    __slots__ = (
        "_exon_indices",
        "_intron_indices",
        "_5utr_indices",
        "_3utr_indices",
        "stop_codons",
        "all_stop_codons",
        "start_codons",
        "all_start_codons",
        "five_prime_utr",
        "three_prime_utr",
        "gene_start",
        "gene_end",
        "strand",
        "gene_type",
        "utr_indices",
        "chrom",
        "genome",
        "_full_mrna",
    )

    def __init__(self, seq=None):
        self._exon_indices = []
        self._intron_indices = []
        self._5utr_indices = []
        self._3utr_indices = []
        self.stop_codons = []
        self.all_stop_codons = []
        self.start_codons = []
        self.all_start_codons = []
        self.five_prime_utr = ""
        self.three_prime_utr = ""
        self.gene_start = None
        self.gene_end = None
        self.strand = None
        self.gene_type = None
        self.utr_indices = []
        self.chrom = None
        self.genome = None
        self._full_mrna = None

        if seq is not None:
            self._exon_indices = [(0, len(seq))]
            self.gene_start = 0
            self.gene_end = len(seq)
            self.strand = StrandType.POS
            self.gene_type = GeneType.UNKNOWN
            self._full_mrna = seq

    @property
    def full_mrna(self):
        if (
            self._full_mrna is None
            and self.chrom is not None
            and self.gene_start is not None
            and self.gene_end is not None
        ):
            # Lazy import to avoid a circular import at module load time.
            from .read_human_genome import fetch_full_mrna

            self._full_mrna = fetch_full_mrna(self.genome, self.chrom, self.gene_start, self.gene_end, self.strand)
        return self._full_mrna

    @full_mrna.setter
    def full_mrna(self, value):
        self._full_mrna = value

    def _get_sequence_slice(self, start, end):
        full_mrna = self.full_mrna
        if not full_mrna or self.gene_start is None or self.gene_end is None:
            return ""
        if self.strand == StrandType.POS:
            rel_start = start - self.gene_start
            rel_end = end - self.gene_start
        else:
            rel_start = self.gene_end - end
            rel_end = self.gene_end - start
        return full_mrna[rel_start:rel_end]

    def add_exon_indices(self, start, end):
        self._exon_indices.append((start, end))

    def add_intron_indices(self, start, end):
        self._intron_indices.append((start, end))

    def add_5utr_indices(self, start, end):
        self._5utr_indices.append((start, end))

    def add_3utr_indices(self, start, end):
        self._3utr_indices.append((start, end))

    @property
    def exons(self):
        return [self._get_sequence_slice(s, e) for s, e in self._exon_indices]

    @property
    def introns(self):
        return [self._get_sequence_slice(s, e) for s, e in self._intron_indices]

    @property
    def sense_utr5(self):
        return [self._get_sequence_slice(s, e) for s, e in self._5utr_indices]

    @property
    def sense_utr3(self):
        return [self._get_sequence_slice(s, e) for s, e in self._3utr_indices]

    def __repr__(self):
        mrna_len = len(self._full_mrna) if self._full_mrna else 0
        loaded = "loaded" if self._full_mrna is not None else "lazy"
        strand_sym = "+" if self.strand == StrandType.POS else "-" if self.strand == StrandType.NEG else "?"
        gene_type_name = self.gene_type.name if self.gene_type is not None else None
        return (
            f"LazyLocusInfo(strand={strand_sym}, type={gene_type_name}, "
            f"coords={self.gene_start}-{self.gene_end}, mrna={loaded}({mrna_len}), "
            f"exons={len(self._exon_indices)}, introns={len(self._intron_indices)})"
        )
