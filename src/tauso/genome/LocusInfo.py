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


def _coding_exon_intervals(exon_intervals, utr_intervals):
    """Coding-exon intervals = exon intervals with any UTR overlap removed.

    Inputs and output are 0-based half-open genomic [start, end) intervals; the result
    is sorted by ascending genomic start. This is purely genomic interval arithmetic and
    so is strand-agnostic. It is the genomic-coordinate counterpart of the structure
    pipeline's ``is_exon & ~is_utr`` test (a position lies in a returned interval iff it
    is exonic and not in a UTR); for a protein-coding transcript that span is exactly the
    coding sequence, from the start codon through the stop codon.
    """
    utr_intervals = sorted(utr_intervals)
    coding_intervals = []
    for exon_start, exon_end in sorted(exon_intervals):
        fragments = [(exon_start, exon_end)]
        for utr_start, utr_end in utr_intervals:
            remaining = []
            for frag_start, frag_end in fragments:
                if utr_end <= frag_start or utr_start >= frag_end:  # disjoint
                    remaining.append((frag_start, frag_end))
                    continue
                if frag_start < utr_start:
                    remaining.append((frag_start, utr_start))
                if utr_end < frag_end:
                    remaining.append((utr_end, frag_end))
            fragments = remaining
        coding_intervals.extend(frag for frag in fragments if frag[1] > frag[0])
    return sorted(coding_intervals)


class CanonicalCds:
    """A locus's canonical coding sequence and the pre-mRNA -> CDS position map.

    Built by :meth:`from_locus` from the same exon/UTR intervals that define the
    ``sense_cds`` flag, so codon features and that flag share one transcript model.

    ``sequence`` is the spliced CDS (start codon through stop codon) in 5'->3'
    orientation. ``premrna_intervals`` is a list of (premrna_start, premrna_end,
    cds_offset) coding spans in pre-mRNA coordinates (an index into the gene-span
    ``full_mrna``), in biological 5'->3' order, where cds_offset is the cumulative
    CDS length preceding the span. Both are empty for non-coding genes.
    """

    def __init__(self, sequence, premrna_intervals):
        self.sequence = sequence
        self.premrna_intervals = premrna_intervals

    def to_cds_index(self, premrna_index):
        """0-based index into ``sequence`` for a pre-mRNA position, or None if the
        position is not inside the CDS. Strand-agnostic: ``premrna_intervals`` are
        already in pre-mRNA coordinates, so the strand was resolved when they were built.
        """
        for premrna_start, premrna_end, cds_offset in self.premrna_intervals:
            if premrna_start <= premrna_index < premrna_end:
                return cds_offset + (premrna_index - premrna_start)
        return None

    @classmethod
    def from_locus(cls, locus):
        """Build the canonical CDS from a locus. Non-coding genes (anything other than a
        PROTEIN_CODING biotype) get an empty CDS, mirroring the ``sense_cds`` biotype guard.
        """
        if locus.gene_type != GeneType.PROTEIN_CODING:
            return cls("", [])

        coding_intervals = _coding_exon_intervals(locus._exon_indices, locus._5utr_indices + locus._3utr_indices)
        # Walk the coding exons in biological 5'->3' order so the sequence reads
        # start->stop and cds_offset accumulates in translation order. The interval
        # arithmetic above is genomic (ascending); only the iteration order is flipped
        # here for the minus strand.
        is_minus_strand = locus.strand == StrandType.NEG
        if is_minus_strand:
            coding_intervals = sorted(coding_intervals, reverse=True)

        gene_start, gene_end = locus.gene_start, locus.gene_end
        premrna_intervals = []
        sequence_parts = []
        cds_offset = 0
        for exon_start, exon_end in coding_intervals:
            # Convert the genomic coding span to pre-mRNA coordinates (an index into
            # full_mrna, which the loader already stored 5'->3'). This is where the
            # strand is resolved.
            if is_minus_strand:
                premrna_start = gene_end - exon_end
                premrna_end = gene_end - exon_start
            else:
                premrna_start = exon_start - gene_start
                premrna_end = exon_end - gene_start
            premrna_intervals.append((premrna_start, premrna_end, cds_offset))
            cds_offset += exon_end - exon_start
            sequence_parts.append(locus._get_sequence_slice(exon_start, exon_end))

        return cls("".join(sequence_parts), premrna_intervals)


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

    def __init__(self, seq=None, cds_start=None, cds_end=None):
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
            if (cds_start is None) != (cds_end is None):
                raise ValueError("pass both cds_start and cds_end (0-based, half-open), or neither")
            if cds_start is not None:
                self._annotate_single_exon_cds(cds_start, cds_end)

    def _annotate_single_exon_cds(self, cds_start, cds_end):
        """Annotate [cds_start, cds_end) (0-based, half-open) as the coding region of a single-exon
        custom transcript -- marks it PROTEIN_CODING with flanking 5'/3'UTRs and start/stop codons,
        so struct_sense_in_5utr/cds/3utr and the start/stop-distance features come out correct."""
        n = len(self.full_mrna)
        if not (0 <= cds_start < cds_end <= n) or (cds_end - cds_start) % 3 != 0:
            raise ValueError(f"cds {(cds_start, cds_end)} is not a valid in-frame coding span for length {n}")
        self.gene_type = GeneType.PROTEIN_CODING
        self._5utr_indices = [(0, cds_start)] if cds_start > 0 else []
        self._3utr_indices = [(cds_end, n)] if cds_end < n else []
        self.utr_indices = self._5utr_indices + self._3utr_indices
        self.five_prime_utr = self.full_mrna[:cds_start]
        self.three_prime_utr = self.full_mrna[cds_end:]
        self.start_codons = self.all_start_codons = [(cds_start, cds_start + 3)]
        self.stop_codons = self.all_stop_codons = [(cds_end - 3, cds_end)]

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

    @property
    def cds_sequence(self):
        """Spliced canonical CDS (start codon .. stop codon), 5'->3'; '' if non-coding."""
        return CanonicalCds.from_locus(self).sequence

    @property
    def cds_premrna_intervals(self):
        """Coding spans as (premrna_start, premrna_end, cds_offset); see CanonicalCds."""
        return CanonicalCds.from_locus(self).premrna_intervals

    def premrna_to_cds_index(self, premrna_index):
        """0-based CDS index for a pre-mRNA position, or None if outside the CDS."""
        return CanonicalCds.from_locus(self).to_cds_index(premrna_index)

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

    @property
    def cds_sequence(self):
        """Spliced canonical CDS (start codon .. stop codon), 5'->3'; '' if non-coding."""
        return CanonicalCds.from_locus(self).sequence

    @property
    def cds_premrna_intervals(self):
        """Coding spans as (premrna_start, premrna_end, cds_offset); see CanonicalCds."""
        return CanonicalCds.from_locus(self).premrna_intervals

    def premrna_to_cds_index(self, premrna_index):
        """0-based CDS index for a pre-mRNA position, or None if outside the CDS."""
        return CanonicalCds.from_locus(self).to_cds_index(premrna_index)

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
