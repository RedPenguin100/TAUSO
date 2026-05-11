import logging
from collections import defaultdict

from ..data.data import load_db, load_genome
from ..timer import Timer
from ..util import get_antisense_u
from .LocusInfo import GeneType, LocusInfo, StrandType

logger = logging.getLogger(__name__)


def get_locus_to_data_dict(include_introns=True, gene_subset=None, genome="GRCh38"):
    with Timer() as t:
        db = load_db(genome)
    logger.debug(f"[Get_Locus] Loaded annotation database in: {t.elapsed_time}s")
    with Timer() as t:
        fasta_dict = load_genome(genome)
    logger.debug(f"[Get_Locus] Loaded fasta dict in: {t.elapsed_time}s")

    locus_to_data = defaultdict(LocusInfo)

    basic_features = ("exon", "gene", "UTR")
    feature_types = ("exon", "gene", "UTR", "intron") if include_introns else basic_features

    iterator = []

    # If we have specific genes, don't scan the whole genome.
    target_names = None

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

        if target_names and locus_tag not in target_names:
            continue

        locus_info = locus_to_data[locus_tag] # defaultdict will assign the default if missing

        locus_info.strand = StrandType.from_string(feature.strand)

        if feature.featuretype == "exon":
            locus_info.add_exon_indices(feature.start - 1, feature.end)

        elif feature.featuretype == "intron" and include_introns:
            locus_info.add_intron_indices(feature.start - 1, feature.end)

        elif feature.featuretype == "gene":
            # We ONLY pull the sequence into memory at the Gene level
            seq = str(fasta_dict[chrom][feature.start - 1 : feature.end])
            if locus_info.strand == StrandType.NEG:
                seq = get_antisense_u(seq)
            else:
                seq = seq.upper()

            locus_info.gene_start = feature.start - 1
            locus_info.gene_end = feature.end
            locus_info.full_mrna = seq

            raw_type_list = feature.attributes.get("gene_type", feature.attributes.get("gene_biotype", []))
            raw_type_str = raw_type_list[0] if raw_type_list else "unannotated"

            locus_info.gene_type = GeneType.from_string(raw_type_str)

        elif "UTR" in feature.featuretype:
            locus_info.utr_indices.append((feature.start - 1, feature.end))
        else:
            logger.debug(f"[Get_Locus] Feature type: {feature.featuretype}")

    # Final cleanup: Reverse coordinates for the negative strand so
    # exon_indices go biologically 5' to 3'.
    for locus_tag, locus_info in locus_to_data.items():
        locus_info._exon_indices.sort()
        locus_info.utr_indices.sort()
        if include_introns:
            locus_info._intron_indices.sort()

        if locus_info.strand == StrandType.NEG:
            locus_info._exon_indices.reverse()
            if include_introns:
                locus_info._intron_indices.reverse()

    return locus_to_data
