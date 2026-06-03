import pytest

from tauso.genome.LocusInfo import GeneType, StrandType
from tauso.genome.read_human_genome import get_locus_to_data_dict

# ── Ground truth ───────────────────────────────────────────────────────────────

GROUND_TRUTH = {
    "BRCA1": {
        "strand": StrandType.NEG,
        "gene_type": GeneType.PROTEIN_CODING,
        "chrom": "chr17",
        "bio_five_prime": 43170245,
        "bio_three_prime": 43044294,
    },
    "BRCA2": {
        "strand": StrandType.POS,
        "gene_type": GeneType.PROTEIN_CODING,
        "chrom": "chr13",
        "bio_five_prime": 32315085,
        "bio_three_prime": 32400268,
    },
    "ACTB": {
        "strand": StrandType.NEG,
        "gene_type": GeneType.PROTEIN_CODING,
        "chrom": "chr7",
        "bio_five_prime": 5563902,
        "bio_three_prime": 5526408,
    },
    "TP53": {
        "strand": StrandType.NEG,
        "gene_type": GeneType.PROTEIN_CODING,
        "chrom": "chr17",
        "bio_five_prime": 7687538,
        "bio_three_prime": 7661778,
    },
    "MALAT1": {
        "strand": StrandType.POS,
        "gene_type": GeneType.LNCRNA,
        "chrom": "chr11",
        "bio_five_prime": 65497687,
        "bio_three_prime": 65506516,
        "canonical_n_exons": 1,
    },
    "NEAT1": {
        "strand": StrandType.POS,
        "gene_type": GeneType.LNCRNA,
        "chrom": "chr11",
        "bio_five_prime": 65422773,
        "bio_three_prime": 65445540,
    },
    "XIST": {
        "strand": StrandType.NEG,
        "gene_type": GeneType.LNCRNA,
        "chrom": "chrX",
        "bio_five_prime": 73852723,
        "bio_three_prime": 73820648,
    },
    "GAPDH": {
        "strand": StrandType.POS,
        "gene_type": GeneType.PROTEIN_CODING,
        "chrom": "chr12",
        "bio_five_prime": 6534511,
        "bio_three_prime": 6538374,
    },
    "RBFOX1": {
        "strand": StrandType.POS,
        "gene_type": GeneType.PROTEIN_CODING,
        "chrom": "chr16",
        "bio_five_prime": 5239801,
        "bio_three_prime": 7713340,
    },
}
REGRESSION_GENES = list(GROUND_TRUTH.keys())
UTR_GENES = ["BRCA1", "BRCA2", "TP53", "GAPDH", "ACTB"]

# ── Fixtures ───────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def locus_data_non_canonical():
    """All transcripts — default mode."""
    return get_locus_to_data_dict(
        include_introns=True,
        gene_subset=REGRESSION_GENES,
        genome="GRCh38",
        canonical_only=False,
    )


@pytest.fixture(scope="module")
def locus_data_canonical():
    """Canonical transcript only — matches literature for single-exon genes."""
    return get_locus_to_data_dict(
        include_introns=True,
        gene_subset=REGRESSION_GENES,
        genome="GRCh38",
        canonical_only=True,
    )


# ── Ground Truth Parameterized Test ────────────────────────────────────────────


@pytest.mark.parametrize("gene, truth", GROUND_TRUTH.items())
def test_ground_truth_properties(locus_data_non_canonical, gene, truth):
    """Combines presence, strand, gene_type, and biological coordinate validations."""
    assert gene in locus_data_non_canonical, f"{gene} missing from result"
    info = locus_data_non_canonical[gene]

    assert info.strand == truth["strand"], f"Expected {truth['strand']}, got {info.strand}"
    assert info.gene_type == truth["gene_type"], f"Expected {truth['gene_type']}, got {info.gene_type}"

    # Verify biological 5' and 3' ends
    first_exon = info._exon_indices[0]
    last_exon = info._exon_indices[-1]

    actual_five_prime = first_exon[0] if info.strand == StrandType.POS else first_exon[1]
    actual_three_prime = last_exon[1] if info.strand == StrandType.POS else last_exon[0]

    assert actual_five_prime == truth["bio_five_prime"], f"bio 5' end {actual_five_prime} != {truth['bio_five_prime']}"
    assert actual_three_prime == truth["bio_three_prime"], (
        f"bio 3' end {actual_three_prime} != {truth['bio_three_prime']}"
    )


def test_chrX_not_filtered(locus_data_canonical):
    """chrX must survive the canonical chrom filter."""
    assert "XIST" in locus_data_canonical, "XIST on chrX was incorrectly filtered out"


def test_no_duplicate_exon_coords(locus_data_non_canonical):
    """After dedup, no two exons should have identical (start, end)."""
    for gene, info in locus_data_non_canonical.items():
        coords = info._exon_indices
        assert len(coords) == len(set(coords)), f"{gene}: duplicate exon coordinates found"


def test_no_overlapping_exons(locus_data_non_canonical):
    """Exons must not overlap after merge and sort."""
    for gene, info in locus_data_non_canonical.items():
        exons = sorted(info._exon_indices)
        for i in range(len(exons) - 1):
            assert exons[i][1] <= exons[i + 1][0], (
                f"{gene}: overlapping exons at index {i}: {exons[i]} vs {exons[i + 1]}"
            )


def test_introns_are_exactly_exon_gaps(locus_data_non_canonical):
    """Derived introns must equal the gaps between merged sorted exons."""
    for gene, info in locus_data_non_canonical.items():
        exons = sorted(info._exon_indices)
        derived = [(exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1) if exons[i][1] < exons[i + 1][0]]
        assert derived == sorted(info._intron_indices), f"{gene}: introns don't match exon gaps"


def test_intron_count_is_exons_minus_one(locus_data_non_canonical):
    """After merging, introns = exons - 1 for any multi-exon gene."""
    for gene, info in locus_data_non_canonical.items():
        n_exons, n_introns = len(info._exon_indices), len(info._intron_indices)
        if n_exons > 1:
            assert n_introns == n_exons - 1, f"{gene}: expected {n_exons - 1} introns, got {n_introns}"
        else:
            assert n_introns == 0, f"{gene}: single-exon gene should have 0 introns, got {n_introns}"


def test_no_zero_length_introns(locus_data_non_canonical):
    for gene, info in locus_data_non_canonical.items():
        for s, e in info._intron_indices:
            assert e > s, f"{gene}: zero-length intron ({s}, {e})"


# ── General Coords & Sequences ─────────────────────────────────────────────────


def test_gene_coordinates_sane(locus_data_non_canonical):
    for gene, info in locus_data_non_canonical.items():
        assert info.gene_start >= 0
        assert info.gene_end > info.gene_start
        assert len(info.full_mrna) == info.gene_end - info.gene_start, f"{gene}: sequence length != gene span"


def test_gene_contains_all_exons(locus_data_non_canonical):
    for gene, info in locus_data_non_canonical.items():
        for s, e in info._exon_indices:
            assert s >= info.gene_start and e <= info.gene_end, f"{gene}: exon out of gene bounds"


def test_exon_order_is_biological(locus_data_non_canonical):
    """Index 0 must be the biological 5' end."""
    for gene, info in locus_data_non_canonical.items():
        if len(info._exon_indices) < 2:
            continue
        if info.strand == StrandType.NEG:
            assert info._exon_indices[0][0] > info._exon_indices[-1][0], f"{gene}: NEG strand not 5'->3'"
        else:
            assert info._exon_indices[0][0] < info._exon_indices[-1][0], f"{gene}: POS strand not 5'->3'"


def test_sequence_validity(locus_data_non_canonical):
    valid = set("ACGTU")
    for gene, info in locus_data_non_canonical.items():
        assert info.full_mrna == info.full_mrna.upper(), f"{gene}: sequence contains lowercase"
        bad = set(info.full_mrna) - valid
        assert not bad, f"{gene}: unexpected bases {bad}"


def test_large_gene_span(locus_data_non_canonical):
    info = locus_data_non_canonical["RBFOX1"]
    span = info.gene_end - info.gene_start
    assert span > 2_000_000, f"RBFOX1 span {span} looks too small"


def test_gene_boundaries_independent_of_canonical(locus_data_non_canonical, locus_data_canonical):
    for gene in REGRESSION_GENES:
        info_all, info_can = locus_data_non_canonical[gene], locus_data_canonical[gene]
        assert info_all.gene_start == info_can.gene_start
        assert info_all.gene_end == info_can.gene_end
        assert info_all.full_mrna == info_can.full_mrna
        assert info_all.strand == info_can.strand


# ── UTR handling tests ─────────────────────────────────────────────────────────


@pytest.mark.parametrize("gene", UTR_GENES)
def test_utr_presence(locus_data_non_canonical, gene):
    """Protein-coding genes must have UTR intervals successfully assigned."""
    info = locus_data_non_canonical[gene]
    assert info.utr_indices, f"{gene}: utr_indices is empty"
    assert info._5utr_indices, f"{gene}: _5utr_indices is empty"
    assert info._3utr_indices, f"{gene}: _3utr_indices is empty"


@pytest.mark.parametrize("gene", UTR_GENES)
def test_typed_utrs_are_subset_of_utr_indices(locus_data_non_canonical, gene):
    info = locus_data_non_canonical[gene]
    utr_set = set(map(tuple, info.utr_indices))

    for interval in info._5utr_indices:
        assert tuple(interval) in utr_set, f"{gene}: 5' UTR interval {interval} missing from utr_indices"
    for interval in info._3utr_indices:
        assert tuple(interval) in utr_set, f"{gene}: 3' UTR interval {interval} missing from utr_indices"


@pytest.mark.parametrize("gene", UTR_GENES)
def test_utr_bounds_and_coverage(locus_data_non_canonical, gene):
    """UTRs must fall within the gene body, be inside exons, and have length > 0."""
    info = locus_data_non_canonical[gene]
    exon_set = sorted(info._exon_indices)

    for s, e in info.utr_indices:
        assert e > s, f"{gene}: zero-length UTR interval ({s}, {e})"
        assert s >= info.gene_start and e <= info.gene_end, f"{gene}: UTR outside gene bounds"

        covered = any(es <= s and e <= ee for es, ee in exon_set)
        assert covered, f"{gene}: UTR ({s}, {e}) not contained within any exon"


@pytest.mark.parametrize("gene", UTR_GENES)
def test_utr_canonical_vs_all_consistent(locus_data_canonical, gene):
    """The fix must work for canonical_only=True."""
    info = locus_data_canonical[gene]
    assert info._5utr_indices, f"{gene} [canonical]: _5utr_indices empty"
    assert info._3utr_indices, f"{gene} [canonical]: _3utr_indices empty"


# ── Canonical CDS API ───────────────────────────────────────────────────────────

CODING_GENES = ["BRCA1", "BRCA2", "ACTB", "TP53", "GAPDH"]
NONCODING_GENES = ["MALAT1", "NEAT1", "XIST"]
_START_CODONS = {"ATG", "AUG"}
_STOP_CODONS = {"TAA", "TAG", "TGA", "UAA", "UAG", "UGA"}


@pytest.mark.parametrize("gene", CODING_GENES)
def test_cds_sequence_is_valid_orf(locus_data_canonical, gene):
    """The canonical CDS is a clean ORF: multiple of 3, starts ATG, ends on a stop."""
    cds = locus_data_canonical[gene].cds_sequence
    assert cds, f"{gene}: empty canonical CDS"
    assert len(cds) % 3 == 0, f"{gene}: CDS length {len(cds)} not a multiple of 3"
    assert cds[:3] in _START_CODONS, f"{gene}: CDS does not start with a start codon ({cds[:3]})"
    assert cds[-3:] in _STOP_CODONS, f"{gene}: CDS does not end on a stop codon ({cds[-3:]})"


@pytest.mark.parametrize("gene", NONCODING_GENES)
def test_noncoding_genes_have_no_cds(locus_data_canonical, gene):
    """lncRNAs (MALAT1, NEAT1, XIST) have no CDS — the biotype guard mirrors sense_cds."""
    info = locus_data_canonical[gene]
    assert info.cds_sequence == "", f"{gene}: non-coding gene produced a CDS"
    assert info.cds_premrna_intervals == []
    assert info.premrna_to_cds_index(0) is None


@pytest.mark.parametrize("gene", CODING_GENES)
def test_cds_premrna_intervals_reconstruct_cds(locus_data_canonical, gene):
    """Concatenating full_mrna over the coding intervals reproduces cds_sequence,
    and the stored offsets are cumulative coding lengths."""
    info = locus_data_canonical[gene]
    mrna = info.full_mrna
    intervals = info.cds_premrna_intervals
    assert "".join(mrna[s:e] for s, e, _ in intervals) == info.cds_sequence
    offset = 0
    for s, e, cds_off in intervals:
        assert cds_off == offset
        offset += e - s
    assert offset == len(info.cds_sequence)


@pytest.mark.parametrize("gene", CODING_GENES)
def test_premrna_to_cds_index_maps_coding_and_rejects_noncoding(locus_data_canonical, gene):
    """Each coding interval's endpoints map to the expected CDS index; a position
    just 5' of the first coding base (5'UTR or intron) is not in the CDS."""
    info = locus_data_canonical[gene]
    intervals = info.cds_premrna_intervals
    for s, e, cds_off in intervals:
        assert info.premrna_to_cds_index(s) == cds_off
        assert info.premrna_to_cds_index(e - 1) == cds_off + (e - 1 - s)
    first_start = intervals[0][0]
    if first_start > 0:
        assert info.premrna_to_cds_index(first_start - 1) is None


@pytest.mark.parametrize("gene", CODING_GENES)
def test_cds_excludes_utr(locus_data_canonical, gene):
    """No coding-CDS pre-mRNA position falls inside a UTR interval (CDS = exon − UTR)."""
    info = locus_data_canonical[gene]
    gs, ge = info.gene_start, info.gene_end
    neg = info.strand == StrandType.NEG
    utr_premrna = []
    for s, e in info._5utr_indices + info._3utr_indices:
        ps, pe = (ge - e, ge - s) if neg else (s - gs, e - gs)
        utr_premrna.append((ps, pe))
    for s, e, _ in info.cds_premrna_intervals:
        for us, ue in utr_premrna:
            assert e <= us or s >= ue, f"{gene}: CDS interval ({s},{e}) overlaps UTR ({us},{ue})"
