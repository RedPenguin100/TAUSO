import bisect
import logging

from Bio.Seq import Seq

from ..data.data import load_db, load_genome
from ..timer import Timer
from .LocusInfo import LocusInfo

logger = logging.getLogger(__name__)


def get_locus_to_data_dict(include_introns=True, gene_subset=None, genome="GRCh38"):
    with Timer() as t:
        db = load_db(genome)
    logger.debug(f"[Get_Locus] Loaded annotation database in: {t.elapsed_time}s")
    with Timer() as t:
        fasta_dict = load_genome(genome)
    logger.debug(f"[Get_Locus] Loaded fasta dict in: {t.elapsed_time}s")

    locus_to_data = dict()
    locus_to_strand = dict()

    basic_features = ["exon", "gene", "stop_codon", "UTR"]

    feature_types = list(basic_features)

    if include_introns:
        feature_types.append("intron")

    iterator = []

    # If we have specific genes, don't scan the whole genome.
    if gene_subset:
        target_names = set(gene_subset)

        # 2. Iterate ONLY genes (Fast, ~60k items vs millions)
        # Note: If you have Gene IDs (ENSG...), you could do db[id] which is O(1).
        # Since you use names (KLKB1), we scan the genes.
        for gene in db.features_of_type("gene"):
            name_list = gene.attributes.get("gene_name", [])
            if name_list and name_list[0] in target_names:
                # Add the gene itself
                iterator.append(gene)

                # Add its children (exons, introns, etc.)
                # db.children() uses the DB index = INSTANT
                children = db.children(gene, featuretype=feature_types, order_by="start")
                iterator.extend(list(children))
    else:
        iterator = db.features_of_type(feature_types, order_by="start")

    for feature in iterator:
        chrom = feature.seqid
        if "chrM" == chrom:
            continue

        locus_tags = feature.attributes["gene_name"]
        if len(locus_tags) != 1:
            raise ValueError(f"Multiple loci: {locus_tags}")
        locus_tag = locus_tags[0]

        if gene_subset is not None and locus_tag not in gene_subset:
            continue

        if locus_tag not in locus_to_data:
            locus_info = LocusInfo()
            locus_to_data[locus_tag] = locus_info
        else:
            locus_info = locus_to_data[locus_tag]

        if feature.featuretype == "exon":
            # Just save the coordinates. The actual sequence is derived lazily.
            bisect.insort(locus_info.exon_indices, (feature.start - 1, feature.end))
            locus_to_strand[locus_tag] = feature.strand

        elif feature.featuretype == "intron" and include_introns:
            bisect.insort(locus_info.intron_indices, (feature.start - 1, feature.end))
            locus_to_strand[locus_tag] = feature.strand

        elif feature.featuretype == "gene":
            # We ONLY pull the sequence into memory at the Gene level
            seq = Seq(str(fasta_dict[chrom][feature.start - 1 : feature.end]))
            if feature.strand == "-":
                seq = seq.reverse_complement()

            locus_info.strand = feature.strand
            locus_info.gene_start = feature.start - 1
            locus_info.gene_end = feature.end
            locus_info.full_mrna = str(seq).upper()
            locus_to_strand[locus_tag] = feature.strand

            raw_type = feature.attributes.get("gene_type", feature.attributes.get("gene_biotype", ["unannotated"]))
            locus_info.gene_type = raw_type[0] if raw_type else "unannotated"

        elif "UTR" in feature.featuretype:
            bisect.insort(locus_info.utr_indices, (feature.start - 1, feature.end))

        elif feature.featuretype == "stop_codon":
            locus_info.stop_codons.append((feature.start, feature.end))

        else:
            logger.debug(f"[Get_Locus] Feature type: {feature.featuretype}")

    # Final cleanup: Reverse coordinates for the negative strand so
    # exon_indices go biologically 5' to 3'.
    for locus_tag, locus_info in locus_to_data.items():
        if locus_to_strand.get(locus_tag) == "-":
            locus_info.exon_indices.reverse()
            if include_introns:
                locus_info.intron_indices.reverse()

    return locus_to_data
