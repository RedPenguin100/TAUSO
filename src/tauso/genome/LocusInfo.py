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


def _subtract_utr_from_exons(exon_intervals, utr_intervals):
    """Coding-exon intervals = exon intervals with any UTR overlap removed.

    All intervals are 0-based half-open genomic [start, end). Returns the coding
    sub-intervals sorted by ascending genomic start. This is the genomic-coordinate
    counterpart of the structure pipeline's ``is_exon & ~is_utr`` test, so a position
    lies in a returned interval iff it is exonic and not in a UTR.
    """
    utrs = sorted(utr_intervals)
    coding = []
    for es, ee in sorted(exon_intervals):
        segments = [(es, ee)]
        for us, ue in utrs:
            nxt = []
            for s, e in segments:
                if ue <= s or us >= e:  # no overlap
                    nxt.append((s, e))
                    continue
                if s < us:
                    nxt.append((s, us))
                if ue < e:
                    nxt.append((ue, e))
            segments = nxt
        coding.extend(seg for seg in segments if seg[1] > seg[0])
    return sorted(coding)


def _coding_exon_intervals_biological(locus):
    """Coding-exon genomic intervals of the canonical transcript, in biological
    5'->3' order. Empty unless the gene is protein-coding (mirrors the biotype
    guard on ``sense_cds``, so lncRNAs such as MALAT1 yield no CDS)."""
    if locus.gene_type != GeneType.PROTEIN_CODING:
        return []
    coding = _subtract_utr_from_exons(locus._exon_indices, locus._5utr_indices + locus._3utr_indices)
    if locus.strand == StrandType.NEG:
        coding = sorted(coding, reverse=True)
    return coding


def compute_cds_sequence(locus):
    """Spliced canonical CDS (ATG..stop), 5'->3'. '' for non-coding genes."""
    return "".join(locus._get_sequence_slice(s, e) for s, e in _coding_exon_intervals_biological(locus))


def compute_cds_premrna_intervals(locus):
    """Coding-exon spans in pre-mRNA index space as (premrna_start, premrna_end,
    cds_offset) tuples in biological order, where pre-mRNA index = position in
    ``full_mrna`` (the gene-span sequence, 5'->3') and cds_offset is the cumulative
    CDS length before this span. Used to map an ASO's pre-mRNA position to a CDS index.
    """
    out = []
    offset = 0
    gs, ge = locus.gene_start, locus.gene_end
    neg = locus.strand == StrandType.NEG
    for s, e in _coding_exon_intervals_biological(locus):
        pstart, pend = (ge - e, ge - s) if neg else (s - gs, e - gs)
        out.append((pstart, pend, offset))
        offset += e - s
    return out


def map_premrna_to_cds_index(premrna_idx, cds_premrna_intervals):
    """0-based index into the CDS for a pre-mRNA position, or None if not in the CDS.

    ``cds_premrna_intervals`` is the output of :func:`compute_cds_premrna_intervals`.
    """
    for pstart, pend, offset in cds_premrna_intervals:
        if pstart <= premrna_idx < pend:
            return offset + (premrna_idx - pstart)
    return None


class _CdsMixin:
    """Canonical-CDS accessors shared by LocusInfo and LazyLocusInfo (both expose
    the same exon/UTR interval attributes and ``_get_sequence_slice``)."""

    __slots__ = ()

    @property
    def cds_sequence(self):
        return compute_cds_sequence(self)

    @property
    def cds_premrna_intervals(self):
        return compute_cds_premrna_intervals(self)

    def premrna_to_cds_index(self, premrna_idx):
        """0-based CDS index for a pre-mRNA position, or None if outside the CDS."""
        return map_premrna_to_cds_index(premrna_idx, compute_cds_premrna_intervals(self))


class LocusInfo(_CdsMixin):
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


class LazyLocusInfo(_CdsMixin):
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
