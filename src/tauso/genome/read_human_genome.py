import logging
from collections import defaultdict

from ..data.data import load_db, load_genome
from ..timer import Timer
from ..util import get_antisense_u
from .LocusInfo import GeneType, LocusInfo, StrandType

logger = logging.getLogger(__name__)

import logging
from collections import defaultdict

from ..data.data import load_db, load_genome
from ..timer import Timer
from ..util import get_antisense_u
from .LocusInfo import GeneType, LocusInfo, StrandType

logger = logging.getLogger(__name__)


def _merge_overlapping_exons(exon_list):
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
    sorted_exons = sorted(exon_list)
    introns = []
    for i in range(len(sorted_exons) - 1):
        intron_start = sorted_exons[i][1]
        intron_end = sorted_exons[i + 1][0]
        if intron_start < intron_end:
            introns.append((intron_start, intron_end))
    return introns


def _canonical_iter(db, canonical_transcripts, child_feature_types):
    for tx_id in canonical_transcripts:
        try:
            tx = db[tx_id]
            yield from db.children(tx, featuretype=child_feature_types, order_by="start")
        except Exception:
            logger.warning(f"[Get_Locus] Could not fetch transcript {tx_id}")


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
    """
    GFF3-compatible version of get_locus_to_data_dict.

    Key differences from the GTF version:
    - Gene name resolution tries GFF3 attribute names first (Name, gene)
    - Gene biotype resolution tries GFF3 attribute names first (biotype, gene_biotype)
    - child_feature_types includes 'CDS' instead of 'UTR' (GFF3 represents UTRs
      as five_prime_UTR / three_prime_UTR, so we cast on the "UTR" substring check below)
    - Canonical transcript tag check is the same ('Ensembl_canonical') but falls
      back gracefully since NCBI GFF3s don't have this tag at all
    - Duplicate name handling: GFF3 genes on alt loci often share a Name with the
      primary locus, so the existing skip-on-duplicate logic is intentionally kept
    """
    CANONICAL_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

    # GFF3 UTR feature type names differ from GTF
    UTR_FEATURE_TYPES = (
        "five_prime_UTR",
        "three_prime_UTR",
        "five_prime_utr",
        "three_prime_utr",  # some sources lowercase
        "UTR",
    )  # fallback / converted files

    # exon + all known UTR spellings — CDS deliberately excluded because
    # we derive introns from exon gaps rather than trusting CDS records
    child_feature_types = ("exon",) + UTR_FEATURE_TYPES

    with Timer() as t:
        db = load_db(genome)
    logger.debug(f"[Get_Locus] Loaded annotation database in: {t.elapsed_time}s")
    with Timer() as t:
        fasta_dict = load_genome(genome)
    logger.debug(f"[Get_Locus] Loaded fasta dict in: {t.elapsed_time}s")

    locus_to_data = defaultdict(LocusInfo)
    target_names = set(gene_subset) if gene_subset else None

    for gene in db.features_of_type("gene", order_by="start"):
        chrom = gene.seqid

        if chrom not in CANONICAL_CHROMS:
            continue

        locus_tag = _get_gene_name_gff(gene)

        if target_names and locus_tag not in target_names:
            continue

        if locus_tag in locus_to_data:
            logger.warning(
                f"[Get_Locus] Duplicate gene {locus_tag} on {chrom} "
                f"(strand={gene.strand}, start={gene.start}), "
                f"already seen on strand={locus_to_data[locus_tag].strand}. Skipping."
            )
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
            for transcript in db.children(gene, featuretype=("mRNA", "transcript"), order_by="start"):
                tags = transcript.attributes.get("tag", [])

                # Check if 'Ensembl_canonical' is in any of the tag strings
                if any("Ensembl_canonical" in t for t in tags):
                    canonical_transcripts.add(transcript.id)

            if not canonical_transcripts:
                logger.warning(f"[Get_Locus] No canonical transcript for {locus_tag}, using all")
                feature_iter = db.children(gene, featuretype=child_feature_types, order_by="start")
            else:
                feature_iter = _canonical_iter(db, canonical_transcripts, child_feature_types)
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
            elif "UTR" in ft or "utr" in ft:
                # catches five_prime_UTR, three_prime_UTR, five_prime_utr, UTR, etc.
                locus_info.utr_indices.append((feature.start - 1, feature.end))
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
    CANONICAL_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

    with Timer() as t:
        db = load_db(genome)
    logger.debug(f"[Get_Locus] Loaded annotation database in: {t.elapsed_time}s")
    with Timer() as t:
        fasta_dict = load_genome(genome)
    logger.debug(f"[Get_Locus] Loaded fasta dict in: {t.elapsed_time}s")

    locus_to_data = defaultdict(LocusInfo)
    target_names = set(gene_subset) if gene_subset else None
    child_feature_types = ("exon", "UTR", "five_prime_UTR", "three_prime_UTR")

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
            logger.warning(
                f"[Get_Locus] Duplicate gene {locus_tag} on {chrom} "
                f"(strand={gene.strand}, start={gene.start}), "
                f"already seen on strand={locus_to_data[locus_tag].strand}. Skipping."
            )
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

            elif ft in ("five_prime_UTR", "five_prime_utr"):
                locus_info.add_five_prime_utr_indices(feature.start - 1, feature.end)
                locus_info.utr_indices.append((feature.start - 1, feature.end))  # keep legacy field

            elif ft in ("three_prime_UTR", "three_prime_utr"):
                locus_info.add_three_prime_utr_indices(feature.start - 1, feature.end)
                locus_info.utr_indices.append((feature.start - 1, feature.end))

            elif "UTR" in ft or "utr" in ft:
                # Untyped UTR — goes into legacy utr_indices only, can't assign to 5'/3'
                logger.debug(f"[Get_Locus] Untyped UTR for {locus_tag}, cannot assign to 5' or 3'")
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


def get_locus_to_data_dict_old(include_introns=True, gene_subset=None, genome="GRCh38"):
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
        iterator = db.features_of_type(feature_types, order_by="featuretype")
    logger.debug("[Get_Locus] Beginning iteration")
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

        locus_info = locus_to_data[locus_tag]  # defaultdict will assign the default if missing

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

            raw_type_list = feature.attributes.get("gene_type", feature.attributes.get("gene_biotype", ()))
            raw_type_str = raw_type_list[0] if raw_type_list else "unannotated"

            locus_info.gene_type = GeneType.from_string(raw_type_str)

        elif "UTR" in feature.featuretype:
            locus_info.utr_indices.append((feature.start - 1, feature.end))
        else:
            logger.debug(f"[Get_Locus] Feature type: {feature.featuretype}")
    logger.debug("[Get_Locus iteration done]")

    # Final cleanup: Reverse coordinates for the negative strand so
    # exon_indices go biologically 5' to 3'.
    for _, locus_info in locus_to_data.items():
        locus_info._exon_indices.sort()
        locus_info.utr_indices.sort()
        locus_info._five_prime_utr_indices.sort()  # NEW
        locus_info._three_prime_utr_indices.sort()  # NEW
        if include_introns:
            locus_info._intron_indices.sort()

        if locus_info.strand == StrandType.NEG:
            locus_info._exon_indices.reverse()
            locus_info._five_prime_utr_indices.reverse()  # NEW
            locus_info._three_prime_utr_indices.reverse()  # NEW
            if include_introns:
                locus_info._intron_indices.reverse()

    return locus_to_data
