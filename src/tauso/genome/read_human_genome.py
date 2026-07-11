import logging
import os
import pickle
from collections import defaultdict

from ..data.data import get_paths, load_genome, load_gff_db, load_gtf_db
from ..timer import Timer
from ..util import get_antisense_u
from .LocusInfo import GeneType, LazyLocusInfo, LocusInfo, StrandType

# Per-process FASTA handle cache, keyed by (genome, pid). pyfaidx handles are not
# safe to share across forked workers (shared fd offset), so each process opens
# its own; keying on pid makes this correct after a fork.
_PROCESS_FASTA: dict = {}


def fetch_full_mrna(genome, chrom, gene_start, gene_end, strand):
    """Fetch a gene's full genomic span (introns included) from the FASTA and
    apply the same strand transform as the eager loader. Used by LazyLocusInfo
    to materialize `full_mrna` on demand."""
    key = (genome, os.getpid())
    fasta = _PROCESS_FASTA.get(key)
    if fasta is None:
        fasta = load_genome(genome)
        _PROCESS_FASTA[key] = fasta
    seq = str(fasta[chrom][gene_start:gene_end])
    return get_antisense_u(seq) if strand == StrandType.NEG else seq.upper()


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


# Bump when the LazyLocusInfo layout or the builder semantics change in a way
# that makes previously-pickled dicts wrong. A mismatch (or any unpickling error)
# falls back to a rebuild, so a stale cache never silently serves bad coordinates.
_LOCUS_CACHE_VERSION = 4


def _locus_cache_path(genome, include_introns, canonical_only):
    fname = f"{genome}.locus.v{_LOCUS_CACHE_VERSION}.introns{int(include_introns)}.canon{int(canonical_only)}.pkl"
    return os.path.join(get_paths(genome)["dir"], fname)


def _load_locus_cache(cache_path, genome):
    """Wishfully load the full {gene: LazyLocusInfo} dict from a pickle.

    Returns None (so the caller rebuilds) if the pickle is missing, older than
    the gff database it was derived from, or fails to unpickle for any reason
    (e.g. a LazyLocusInfo layout change). Never raises.
    """
    try:
        if not os.path.exists(cache_path):
            return None
        gff_db = get_paths(genome)["gff_db"]
        if os.path.exists(gff_db) and os.path.getmtime(cache_path) < os.path.getmtime(gff_db):
            logger.debug("[Get_Locus] Locus cache %s is older than gff db, rebuilding", cache_path)
            return None
        with open(cache_path, "rb") as fh:
            locus_to_data = pickle.load(fh)
        logger.debug("[Get_Locus] Loaded locus dict from cache %s", cache_path)
        return locus_to_data
    except Exception as e:
        logger.warning("[Get_Locus] Could not load locus cache %s (%s); rebuilding", cache_path, e)
        return None


def _save_locus_cache(cache_path, locus_to_data):
    """Atomically write the full locus dict to the pickle cache. Best-effort:
    a write failure is logged, never fatal."""
    try:
        tmp = f"{cache_path}.tmp{os.getpid()}"
        with open(tmp, "wb") as fh:
            pickle.dump(locus_to_data, fh, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp, cache_path)
        logger.debug("[Get_Locus] Saved locus dict to cache %s", cache_path)
    except Exception as e:
        logger.warning("[Get_Locus] Could not save locus cache %s (%s)", cache_path, e)


def build_locus_cache(genome="GRCh38", include_introns=True, canonical_only=True):
    """Build the full locus dict and persist it to the pickle cache.

    Intended for `setup-genome` so downstream runs (and SLURM tasks) load the
    coordinate dict from pickle instead of re-traversing the gff database. The
    full-dict build inside `get_locus_to_data_dict` writes the cache itself; a
    fresh gff db (e.g. `setup-genome --force`) is newer than any old pickle, so
    the stale cache is bypassed and rewritten. Returns the cache path.
    """
    get_locus_to_data_dict(
        include_introns=include_introns, gene_subset=None, genome=genome, canonical_only=canonical_only
    )
    return _locus_cache_path(genome, include_introns, canonical_only)


def get_locus_to_data_dict(include_introns=True, gene_subset=None, genome="GRCh38", canonical_only=True):
    """GFF3-based genome loader. Builds a {gene_name: LocusInfo} dict.

    Duplicate gene names (e.g. Y_RNA, U6, 5S_rRNA — ncRNAs with hundreds of
    genomic copies) are resolved by keeping the locus with the lowest genomic
    start coordinate, because gffutils iterates with order_by="start" globally
    across all canonical chromosomes.  This is deterministic but not
    biologically motivated; such genes are generally not ASO targets.
    """
    # Wishful load: the cache holds the FULL coordinate-only dict. A subset
    # request is served by filtering it in memory (the global lowest-start dedup
    # is independent of the subset, so this matches a fresh subset build).
    cache_path = _locus_cache_path(genome, include_introns, canonical_only)
    cached = _load_locus_cache(cache_path, genome)
    if cached is not None:
        if not gene_subset:
            return cached
        target_names = set(gene_subset)
        subset = defaultdict(LazyLocusInfo)
        subset.update({name: locus for name, locus in cached.items() if name in target_names})
        return subset

    UTR_FEATURE_TYPES = (
        "five_prime_UTR",
        "three_prime_UTR",
        "five_prime_utr",
        "three_prime_utr",
        "UTR",
    )

    # Canonical transcript information fetching.
    child_feature_types = ("exon", "stop_codon", "start_codon") + UTR_FEATURE_TYPES

    logger.debug("Loading database")
    db = load_gff_db(genome)
    logger.debug("[Get_Locus] Loaded annotation database.")

    locus_to_data = defaultdict(LazyLocusInfo)
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

        locus_info.chrom = chrom
        locus_info.genome = genome
        locus_info.gene_start = gene.start - 1
        locus_info.gene_end = gene.end
        # full_mrna is fetched lazily from the FASTA on first access (LazyLocusInfo).
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
        canonical_stops, canonical_starts = [], []

        for feature in feature_iter:
            key = (feature.featuretype, feature.start, feature.end)
            if key in seen:
                duplicates_skipped += 1
                continue
            seen.add(key)

            ft = feature.featuretype
            if ft == "exon":
                exons_collected.append((feature.start - 1, feature.end))

            elif ft == "stop_codon":
                # A single canonical stop codon is expected for coding transcripts; deviations are handled after the loop.
                canonical_stops.append((feature.start - 1, feature.end))

            elif ft == "start_codon":
                canonical_starts.append((feature.start - 1, feature.end))

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

        # All-transcript stop / start codons (regardless of canonical_only) — used to
        # compute distance to the nearest stop / start codon across any isoform.
        all_stops_seen = set()
        for stop_feat in db.children(gene, featuretype="stop_codon", order_by="start"):
            stop_key = (stop_feat.start, stop_feat.end)
            if stop_key in all_stops_seen:
                continue
            all_stops_seen.add(stop_key)
            locus_info.all_stop_codons.append((stop_feat.start - 1, stop_feat.end))

        all_starts_seen = set()
        for start_feat in db.children(gene, featuretype="start_codon", order_by="start"):
            start_key = (start_feat.start, start_feat.end)
            if start_key in all_starts_seen:
                continue
            all_starts_seen.add(start_key)
            locus_info.all_start_codons.append((start_feat.start - 1, start_feat.end))

        # All-transcript information fetching.
        exons_by_transcript = defaultdict(list)
        for exon in db.children(gene, featuretype="exon", order_by="start"):
            exons_by_transcript[exon.attributes.get("Parent", [None])[0]].append((exon.start - 1, exon.end))
        splice_sites = set()
        for transcript_exons in exons_by_transcript.values():
            for junction_start, junction_end in _derive_introns_from_exons(transcript_exons):
                splice_sites.add(junction_start)
                splice_sites.add(junction_end)
        locus_info.all_splice_junctions = sorted(splice_sites)

        is_coding = locus_info.gene_type == GeneType.PROTEIN_CODING
        locus_info.stop_codon = _standardize_canonical_codons_to_single(canonical_stops, is_coding, locus_tag, "stop")
        locus_info.start_codon = _standardize_canonical_codons_to_single(
            canonical_starts, is_coding, locus_tag, "start"
        )

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

    # Only the full dict is a valid cache; a subset build is partial.
    if not gene_subset:
        _save_locus_cache(cache_path, locus_to_data)

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


def _standardize_canonical_codons_to_single(codons, is_coding, gene_name, kind):
    """For a single coding transcript there should be one stop/start codon, and for a non-coding transcript
    there should be 0. This function handles deviations from that standard format. Not to be confused with
    alternative transcripts, where multiple stop codons can be present."""
    if is_coding:
        if not codons:
            logger.warning(
                "[Get_Locus] %s: protein-coding but 0 canonical %s codons (broken CDS annotation)", gene_name, kind
            )
            return None
        if len(codons) > 1:
            logger.warning(
                "[Get_Locus] %s: %d canonical %s codons (expected 1); keeping the first", gene_name, len(codons), kind
            )
        return codons[0]
    if codons:
        logger.warning("[Get_Locus] %s: non-coding gene has a canonical %s codon; dropping it", gene_name, kind)
    return None


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
