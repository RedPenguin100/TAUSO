from typing import Dict, Iterable, Optional

import gffutils

from ..util import _to_str_seq
from .LocusInfo import GeneType


class GeneCoordinateMapper:
    def __init__(self, db_path: str):
        self.db = gffutils.FeatureDB(db_path)
        self.gene_name_map = {}
        self._cache = {}

        for gene in self.db.features_of_type("gene"):
            try:
                self.gene_name_map[gene["gene_name"][0]] = gene.id
            except KeyError:
                pass

    def _get_gene_and_cds(self, gene_name: str):
        """
        Retrieves the gene and the CDS features (not generic exons)
        of the longest protein-coding transcript.
        """
        if gene_name in self._cache:
            return self._cache[gene_name]

        gene_id = self.gene_name_map.get(gene_name)
        if not gene_id:
            return None, None

        gene = self.db[gene_id]
        transcripts = list(self.db.children(gene, featuretype="transcript"))
        if not transcripts:
            self._cache[gene_name] = (gene, [])
            return gene, []

        # FIX 1: Select Canonical Transcript based on CDS length (Protein Coding potential)
        # We explicitly filter for 'CDS' features to ignore UTR lengths
        def get_cds_len(t):
            return sum(len(c) for c in self.db.children(t, featuretype="CDS"))

        # Filter for transcripts that actually have a CDS
        coding_transcripts = [t for t in transcripts if get_cds_len(t) > 0]

        if not coding_transcripts:
            # Fallback: Gene might be non-coding (lncRNA). Return nothing for CAI purposes.
            self._cache[gene_name] = (gene, [])
            return gene, []

        best_t = max(coding_transcripts, key=get_cds_len)

        # FIX 2: Fetch 'CDS' features, not 'exon'. This excludes UTRs.
        cds_parts = list(self.db.children(best_t, featuretype="CDS", order_by="start"))

        self._cache[gene_name] = (gene, cds_parts)
        return gene, cds_parts

    def get_spliced_index(self, gene_name: str, premrna_idx: int) -> Optional[int]:
        """
        Returns index relative to the Start Codon (ATG).
        Returns None if the ASO is in a UTR or Intron.
        """
        gene, cds_parts = self._get_gene_and_cds(gene_name)
        if not cds_parts:
            return None

        # Determine genomic position of the ASO
        if gene.strand == "+":
            genomic_pos = gene.start + premrna_idx
        else:
            genomic_pos = gene.end - premrna_idx

        spliced_len = 0
        for part in cds_parts:
            # Check if position falls within this CDS block
            if part.start <= genomic_pos <= part.end:
                # Calculate offset inside this block
                if gene.strand == "+":
                    offset = genomic_pos - part.start
                else:
                    offset = part.end - genomic_pos

                # Return cumulative length (previous blocks) + offset
                return spliced_len + offset

            spliced_len += len(part)

        return None  # Position was not in any CDS block (likely Intron or UTR)

    def get_spliced_sequence(self, gene_name: str, pre_mrna_seq: str) -> str:
        """
        Returns the pure Coding Sequence (ATG -> Stop).
        """
        gene, cds_parts = self._get_gene_and_cds(gene_name)
        if not cds_parts:
            return ""

        parts = []
        gene_start = gene.start

        # Since we use 'CDS', this loop builds the sequence strictly from Start Codon
        for part in cds_parts:
            start = part.start - gene_start
            length = part.end - part.start + 1

            # Slice safely from the pre-mRNA sequence string
            if start < len(pre_mrna_seq):
                segment = pre_mrna_seq[max(0, start) : min(len(pre_mrna_seq), start + length)]
                parts.append(segment)

        return "".join(parts)

    def get_gene_coords(self, gene_name: str):
        """
        Returns {chrom, gene_start, gene_end, strand} for the GENE feature.
        """
        gene_id = self.gene_name_map.get(gene_name)
        if not gene_id:
            return None

        # Fetch the gene feature directly
        gene = self.db[gene_id]

        # Convert 1-based GFF to 0-based coordinates
        coords = {
            "chrom": gene.chrom,
            "gene_start": gene.start - 1,
            "gene_end": gene.end,
            "strand": gene.strand,
        }

        return coords


class GeneSequenceRegistry(dict):
    """dict subclass whose repr omits sequence data so pytest failure output stays readable."""

    def __repr__(self) -> str:
        genes = list(self.keys())
        preview = genes[:5]
        suffix = f", ...+{len(genes) - 5} more" if len(genes) > 5 else ""
        return f"GeneSequenceRegistry({len(self)} genes: {preview}{suffix})"


def build_gene_sequence_registry(
    genes: Iterable[str], gene_to_data: Dict, mapper: GeneCoordinateMapper = None
) -> GeneSequenceRegistry:
    """
    Creates a dictionary of
    {gene_name: {'pre_mrna_sequence', 'cds_sequence', 'cds_premrna_intervals'}}
    for fast O(1) lookup during processing.

    The CDS is the canonical-transcript coding sequence taken from the same
    LocusInfo used to compute the sense_cds / sense_exon flags, so codon features
    and those flags share one transcript model. ``cds_premrna_intervals`` maps an
    ASO's pre-mRNA position to a CDS index (see LocusInfo.cds_premrna_intervals).
    ``mapper`` is unused; it is kept for call-site compatibility.
    """
    registry = GeneSequenceRegistry()
    for gene in genes:
        locus = gene_to_data.get(gene)
        if not locus:
            continue

        pre_mrna = _to_str_seq(locus.full_mrna)
        cds = locus.cds_sequence
        cds_intervals = locus.cds_premrna_intervals

        # No canonical CDS: either a non-coding RNA (keep empty) or a custom
        # user-supplied single sequence (treat the whole sequence as the CDS).
        if not cds:
            gene_type = getattr(locus, "gene_type", GeneType.UNKNOWN)
            if gene_type in (GeneType.LNCRNA, GeneType.NON_CODING):
                cds = None
                cds_intervals = []
            else:
                cds = pre_mrna
                cds_intervals = [(0, len(pre_mrna), 0)]

        registry[gene] = {
            "pre_mrna_sequence": pre_mrna,
            "cds_sequence": cds,
            "cds_premrna_intervals": cds_intervals,
        }
    return registry
