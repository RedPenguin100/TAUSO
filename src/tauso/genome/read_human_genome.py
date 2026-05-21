import logging
from collections import defaultdict

from ..data.data import load_genome, load_gff_db, load_gtf_db
from ..timer import Timer
from ..util import get_antisense_u
from .LocusInfo import GeneType, LocusInfo, StrandType

logger = logging.getLogger(__name__)

CANONICAL_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def _get_gene_name_gff(gene):
    """
    GFF3 files are inconsistent about which attribute holds the human-readable
    gene symbol.  Try the most common ones in priority order.

    Ensembl GFF3  → 'Name'
    NCBI GFF3     → 'gene'  (or 'Name')
    Some sources  → 'gene_name' (same as GTF)
    Fallback      → the feature ID so we never silently drop a gene
    """
    for attr in ("Name", "gene", "gene_name"):
        val = gene.attributes.get(attr, [])
        if val:
            return val[0]
    # Last resort: use the ID (e.g. "gene:ENSG00000139618")
    return gene.id


def _get_gene_type_gff(gene):
    """
    GFF3 biotype attribute names vary by source.
    Ensembl GFF3  → 'biotype'
    NCBI GFF3     → 'gene_biotype'
    GTF-style     → 'gene_type' (sometimes preserved in converted GFF3s)
    """
    for attr in ("biotype", "gene_biotype", "gene_type"):
        val = gene.attributes.get(attr, [])
        if val:
            return val[0]
    return "unannotated"


def get_locus_to_data_dict(include_introns=True, gene_subset=None, genome="GRCh38", canonical_only=True):
    """GFF3-based genome loader. Builds a {gene_name: LocusInfo} dict.

    Duplicate gene names (e.g. Y_RNA, U6, 5S_rRNA — ncRNAs with hundreds of
    genomic copies) are resolved by keeping the locus with the lowest genomic
    start coordinate, because gffutils iterates with order_by="start" globally
    across all canonical chromosomes.  This is deterministic but not
    biologically motivated; such genes are generally not ASO targets.
    """
    UTR_FEATURE_TYPES = (
        "five_prime_UTR",
        "three_prime_UTR",
        "five_prime_utr",
        "three_prime_utr",
        "UTR",
    )

    child_feature_types = ("exon",) + UTR_FEATURE_TYPES

    logger.debug("Loading database")
    db = load_gff_db(genome)
    logger.debug("[Get_Locus] Loaded annotation database.")

    logger.debug("[Get_Locus] Loading fasta dict")
    fasta_dict = load_genome(genome)
    logger.debug("[Get_Locus] Loaded fasta dict.")

    locus_to_data = defaultdict(LocusInfo)
    target_names = set(gene_subset) if gene_subset else None
    duplicate_counts: dict[str, int] = {}

    for gene in db.features_of_type("gene", order_by="start"):
        chrom = gene.seqid

        if chrom not in CANONICAL_CHROMS:
            continue

        locus_tag = _get_gene_name_gff(gene)

        if target_names and locus_tag not in target_names:
            continue

        if locus_tag in locus_to_data:
            duplicate_counts[locus_tag] = duplicate_counts.get(locus_tag, 0) + 1
            continue

        locus_info = locus_to_data[locus_tag]
        locus_info.strand = StrandType.from_string(gene.strand)

        seq = str(fasta_dict[chrom][gene.start - 1 : gene.end])
        seq = get_antisense_u(seq) if locus_info.strand == StrandType.NEG else seq.upper()

        locus_info.gene_start = gene.start - 1
        locus_info.gene_end = gene.end
        locus_info.full_mrna = seq
        locus_info.gene_type = GeneType.from_string(_get_gene_type_gff(gene))

        if canonical_only:
            canonical_transcripts = set()
            for transcript in db.children(gene, featuretype=("mRNA", "transcript"), level=1, order_by="start"):
                tags = transcript.attributes.get("tag", [])
                if any("Ensembl_canonical" in t for t in tags):
                    canonical_transcripts.add(transcript.id)

            if not canonical_transcripts:
                logger.warning(f"[Get_Locus] No canonical transcript for {locus_tag}, using all")
                feature_iter = db.children(gene, featuretype=child_feature_types, order_by="start")
            else:
                feature_iter = _canonical_iter(db, canonical_transcripts, child_feature_types)
        else:
            all_transcript_ids = set(tx.id for tx in db.children(gene, featuretype=("mRNA", "transcript"), level=1))
            if all_transcript_ids:
                feature_iter = _canonical_iter(db, all_transcript_ids, child_feature_types)
            else:
                feature_iter = db.children(gene, featuretype=child_feature_types, order_by="start")

        seen = set()
        duplicates_skipped = 0
        exons_collected = []

        for feature in feature_iter:
            key = (feature.featuretype, feature.start, feature.end)
            if key in seen:
                duplicates_skipped += 1
                continue
            seen.add(key)

            ft = feature.featuretype
            if ft == "exon":
                exons_collected.append((feature.start - 1, feature.end))

            elif ft in ("five_prime_UTR", "five_prime_utr"):
                locus_info.add_5utr_indices(feature.start - 1, feature.end)
                locus_info.utr_indices.append((feature.start - 1, feature.end))

            elif ft in ("three_prime_UTR", "three_prime_utr"):
                locus_info.add_3utr_indices(feature.start - 1, feature.end)
                locus_info.utr_indices.append((feature.start - 1, feature.end))

            elif "UTR" in ft or "utr" in ft:
                raise ValueError(
                    f"[Get_Locus] Unrecognized UTR feature type '{ft}' for gene '{locus_tag}' "
                    f"at ({feature.start - 1}, {feature.end}). "
                    f"Add '{ft}' to the five_prime or three_prime dispatch above."
                )

            else:
                logger.debug(f"[Get_Locus] Unexpected feature type: {ft}")

        if duplicates_skipped:
            logger.debug(f"[Get_Locus] {locus_tag}: skipped {duplicates_skipped} duplicate features")

        exons_collected = _merge_overlapping_exons(exons_collected)
        locus_info._exon_indices = []
        for s, e in exons_collected:
            locus_info.add_exon_indices(s, e)

        if include_introns:
            for intron_start, intron_end in _derive_introns_from_exons(exons_collected):
                locus_info.add_intron_indices(intron_start, intron_end)

    if duplicate_counts:
        n_total = sum(duplicate_counts.values())
        examples = ", ".join(list(duplicate_counts)[:5])
        logger.warning(
            f"[Get_Locus] Skipped {n_total} duplicate loci across {len(duplicate_counts)} gene name(s) "
            f"(e.g. {examples}). Kept the locus with the lowest genomic start coordinate."
        )

    logger.debug("[Get_Locus] Iteration done, sorting indices")

    for _, locus_info in locus_to_data.items():
        locus_info._exon_indices.sort()
        locus_info.utr_indices.sort()
        locus_info._5utr_indices.sort()
        locus_info._3utr_indices.sort()
        if include_introns:
            locus_info._intron_indices.sort()

        if locus_info.strand == StrandType.NEG:
            locus_info._exon_indices.reverse()
            locus_info._5utr_indices.reverse()
            locus_info._3utr_indices.reverse()
            if include_introns:
                locus_info._intron_indices.reverse()

    return locus_to_data


def _merge_overlapping_exons(exon_list):
    """Merge overlapping/nested exons into non-overlapping intervals."""
    if not exon_list:
        return []
    sorted_exons = sorted(exon_list)
    merged = [sorted_exons[0]]
    for start, end in sorted_exons[1:]:
        prev_start, prev_end = merged[-1]
        if start <= prev_end:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))
    return merged


def _derive_introns_from_exons(exon_list):
    """Derive introns as gaps between sorted exons."""
    sorted_exons = sorted(exon_list)
    introns = []
    for i in range(len(sorted_exons) - 1):
        intron_start = sorted_exons[i][1]
        intron_end = sorted_exons[i + 1][0]
        if intron_start < intron_end:
            introns.append((intron_start, intron_end))
    return introns


def _canonical_iter(db, canonical_transcripts, child_feature_types):
    """Yield children of canonical transcripts only."""
    for tx_id in canonical_transcripts:
        try:
            tx = db[tx_id]
            yield from db.children(tx, featuretype=child_feature_types, order_by="start")
        except Exception:
            logger.warning(f"[Get_Locus] Could not fetch transcript {tx_id}")


def get_locus_to_data_dict_gtf(include_introns=True, gene_subset=None, genome="GRCh38", canonical_only=False):

    with Timer() as t:
        db = load_gtf_db(genome)
    logger.debug(f"[Get_Locus] Loaded annotation database in: {t.elapsed_time}s")
    with Timer() as t:
        fasta_dict = load_genome(genome)
    logger.debug(f"[Get_Locus] Loaded fasta dict in: {t.elapsed_time}s")

    locus_to_data = defaultdict(LocusInfo)
    target_names = set(gene_subset) if gene_subset else None
    child_feature_types = ("exon", "UTR")
    duplicate_counts: dict[str, int] = {}

    for gene in db.features_of_type("gene", order_by="start"):
        chrom = gene.seqid

        # FIX 1 — scaffold/patch/MT filter
        if chrom not in CANONICAL_CHROMS:
            continue

        name_list = gene.attributes.get("gene_name", [])
        if not name_list:
            continue
        if len(name_list) != 1:
            raise ValueError(f"Multiple loci: {name_list}")
        locus_tag = name_list[0]

        if target_names and locus_tag not in target_names:
            continue

        # FIX 2 — duplicate gene guard (pseudogenes, alt loci, readthroughs)
        if locus_tag in locus_to_data:
            duplicate_counts[locus_tag] = duplicate_counts.get(locus_tag, 0) + 1
            continue

        locus_info = locus_to_data[locus_tag]
        locus_info.strand = StrandType.from_string(gene.strand)

        seq = str(fasta_dict[chrom][gene.start - 1 : gene.end])
        seq = get_antisense_u(seq) if locus_info.strand == StrandType.NEG else seq.upper()

        locus_info.gene_start = gene.start - 1
        locus_info.gene_end = gene.end
        locus_info.full_mrna = seq

        raw_type_list = gene.attributes.get("gene_type", gene.attributes.get("gene_biotype", ()))
        locus_info.gene_type = GeneType.from_string(raw_type_list[0] if raw_type_list else "unannotated")

        # FIX 3 — canonical transcript filter (opt-in)
        # Fetches directly from transcript, bypassing broken Parent attribute on exons
        if canonical_only:
            canonical_transcripts = set()
            for transcript in db.children(gene, featuretype="transcript", order_by="start"):
                tags = transcript.attributes.get("tag", [])
                if "Ensembl_canonical" in tags:
                    canonical_transcripts.add(transcript.id)

            if not canonical_transcripts:
                logger.warning(f"[Get_Locus] No canonical transcript for {locus_tag}, using all")
                feature_iter = db.children(gene, featuretype=child_feature_types, order_by="start")
            else:
                feature_iter = _canonical_iter(db, canonical_transcripts, child_feature_types)
        else:
            feature_iter = db.children(gene, featuretype=child_feature_types, order_by="start")

        # FIX 4 — dedup by exact (type, start, end)
        seen = set()
        duplicates_skipped = 0
        exons_collected = []

        for feature in feature_iter:
            key = (feature.featuretype, feature.start, feature.end)
            if key in seen:
                duplicates_skipped += 1
                continue
            seen.add(key)

            ft = feature.featuretype
            if ft == "exon":
                exons_collected.append((feature.start - 1, feature.end))

            elif ft == "UTR":
                locus_info.utr_indices.append((feature.start - 1, feature.end))

            else:
                logger.debug(f"[Get_Locus] Unexpected feature type: {ft}")

        if duplicates_skipped:
            logger.debug(f"[Get_Locus] {locus_tag}: skipped {duplicates_skipped} duplicate features")

        # FIX 5 — merge overlapping exons FIRST, then derive introns from clean set
        exons_collected = _merge_overlapping_exons(exons_collected)

        locus_info._exon_indices = []
        for s, e in exons_collected:
            locus_info.add_exon_indices(s, e)

        # FIX 6 — derive introns from merged exon gaps, never from db
        if include_introns:
            for intron_start, intron_end in _derive_introns_from_exons(exons_collected):
                locus_info.add_intron_indices(intron_start, intron_end)

    if duplicate_counts:
        n_total = sum(duplicate_counts.values())
        examples = ", ".join(list(duplicate_counts)[:5])
        logger.warning(
            f"[Get_Locus] Skipped {n_total} duplicate loci across {len(duplicate_counts)} gene name(s) "
            f"(e.g. {examples}). Kept the locus with the lowest genomic start coordinate."
        )

    logger.debug("[Get_Locus] Iteration done, sorting indices")

    for _, locus_info in locus_to_data.items():
        locus_info._exon_indices.sort()
        locus_info.utr_indices.sort()
        if include_introns:
            locus_info._intron_indices.sort()

        if locus_info.strand == StrandType.NEG:
            locus_info._exon_indices.reverse()
            if include_introns:
                locus_info._intron_indices.reverse()

    return locus_to_data
