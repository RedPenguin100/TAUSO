from typing import Dict, Iterable

from ..util import _to_str_seq
from .LocusInfo import CanonicalCds, GeneType


class GeneSequenceRegistry(dict):
    """dict subclass whose repr omits sequence data so pytest failure output stays readable."""

    def __repr__(self) -> str:
        genes = list(self.keys())
        preview = genes[:5]
        suffix = f", ...+{len(genes) - 5} more" if len(genes) > 5 else ""
        return f"GeneSequenceRegistry({len(self)} genes: {preview}{suffix})"


def build_gene_sequence_registry(genes: Iterable[str], gene_to_data: Dict) -> GeneSequenceRegistry:
    """
    Creates a dictionary of {gene_name: {'pre_mrna_sequence', 'cds'}} for fast O(1)
    lookup during processing, where 'cds' is a :class:`CanonicalCds`.

    The CDS is the canonical-transcript coding sequence taken from the same LocusInfo
    used to compute the sense_cds / sense_exon flags, so codon features and those flags
    share one transcript model.
    """
    registry = GeneSequenceRegistry()
    for gene in genes:
        locus = gene_to_data.get(gene)
        if not locus:
            continue

        pre_mrna = _to_str_seq(locus.full_mrna)
        cds = CanonicalCds.from_locus(locus)

        # No canonical CDS: a non-coding RNA keeps the empty CDS; a custom user-supplied
        # single sequence (no annotation) is treated as one contiguous CDS.
        if not cds.sequence:
            gene_type = getattr(locus, "gene_type", GeneType.UNKNOWN)
            if gene_type not in (GeneType.LNCRNA, GeneType.NON_CODING):
                cds = CanonicalCds(pre_mrna, [(0, len(pre_mrna), 0)])

        registry[gene] = {"pre_mrna_sequence": pre_mrna, "cds": cds}
    return registry
