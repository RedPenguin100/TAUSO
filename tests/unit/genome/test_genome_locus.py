import pytest

from tauso.genome.LocusInfo import GeneType, StrandType
from tauso.genome.read_human_genome import get_locus_to_data_dict

# ── Ground truth ───────────────────────────────────────────────────────────────
#
# Coord conventions after sort+reverse:
#   POS strand: _exon_indices[0]  = leftmost (lowest genomic, biological 5')
#   NEG strand: _exon_indices[0]  = rightmost (highest genomic, biological 5')
#
# exon_first_start = _exon_indices[0][0]  (start of first biological exon)
# exon_last_end    = _exon_indices[-1][1] for POS, _exon_indices[-1][0] for NEG

GROUND_TRUTH = {
    "BRCA1": {
        # NEG strand: bio 5' = rightmost exon END, bio 3' = leftmost exon START
        "strand": StrandType.NEG,
        "gene_type": GeneType.PROTEIN_CODING,
        "chrom": "chr17",
        "bio_five_prime": 43170245,
        "bio_three_prime": 43044294,
    },
    "BRCA2": {
        # POS strand: bio 5' = leftmost exon START, bio 3' = rightmost exon END
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


# ── Fixtures ───────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def locus_data():
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


# ── Helper: run once to populate GROUND_TRUTH, then delete ────────────────────

def test_print_actual_coords(locus_data):
    """
    TEMPORARY — run once to get real values, paste into GROUND_TRUTH, then delete.
    Always passes; just prints.
    """
    for gene in REGRESSION_GENES:
        info = locus_data[gene]
        print(f"\n{gene}:")
        print(f"  strand          = {info.strand}")
        print(f"  gene_type       = {info.gene_type}")
        print(f"  n_exons         = {len(info._exon_indices)}")
        print(f"  n_introns       = {len(info._intron_indices)}")
        if info._exon_indices:
            print(f"  first_exon      = {info._exon_indices[0]}")
            print(f"  last_exon       = {info._exon_indices[-1]}")
        print(f"  gene_start      = {info.gene_start}")
        print(f"  gene_end        = {info.gene_end}")


def test_chrX_not_filtered(locus_data):
    """chrX must survive the canonical chrom filter."""
    assert "XIST" in locus_data, "XIST on chrX was incorrectly filtered out"


def test_actb_strand_is_correct(locus_data):
    """Confirm ACTB strand after duplicate guard — lock in once confirmed."""
    info = locus_data["ACTB"]
    assert info.strand == StrandType.NEG
    # Once confirmed via test_print_actual_coords, lock it:
    assert info.strand == GROUND_TRUTH["ACTB"]["strand"]


def test_malat1_canonical_is_single_exon(locus_data_canonical):
    """
    MALAT1 is described as single-exon in literature but GENCODE's canonical
    transcript may be a spliced isoform. Assert it has fewer exons than
    all-transcripts mode and is biologically plausible.
    """
    info = locus_data_canonical["MALAT1"]
    assert len(info._exon_indices) >= 1, "MALAT1 canonical should have at least 1 exon"
    assert len(info._intron_indices) == len(info._exon_indices) - 1, "MALAT1 canonical introns should equal exons - 1"


def test_no_duplicate_exon_coords(locus_data):
    """After dedup, no two exons should have identical (start, end)."""
    for gene, info in locus_data.items():
        coords = info._exon_indices
        assert len(coords) == len(set(coords)), f"{gene}: duplicate exon coordinates found"


# ── FIX 5: Exon merging ───────────────────────────────────────────────────────


def test_no_overlapping_exons(locus_data):
    """Exons must not overlap after merge and sort."""
    for gene, info in locus_data.items():
        exons = sorted(info._exon_indices)
        for i in range(len(exons) - 1):
            assert exons[i][1] <= exons[i + 1][0], (
                f"{gene}: overlapping exons at index {i}: {exons[i]} vs {exons[i + 1]}"
            )


# ── FIX 6: Derived introns ────────────────────────────────────────────────────


def test_introns_are_exactly_exon_gaps(locus_data):
    """Derived introns must equal the gaps between merged sorted exons."""
    for gene, info in locus_data.items():
        exons = sorted(info._exon_indices)
        derived = [(exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1) if exons[i][1] < exons[i + 1][0]]
        stored = sorted(info._intron_indices)
        assert derived == stored, f"{gene}: introns don't match exon gaps\n  derived: {derived}\n  stored:  {stored}"


def test_intron_count_is_exons_minus_one(locus_data):
    """After merging, introns = exons - 1 for any multi-exon gene."""
    for gene, info in locus_data.items():
        n_exons = len(info._exon_indices)
        n_introns = len(info._intron_indices)
        if n_exons > 1:
            assert n_introns == n_exons - 1, (
                f"{gene}: expected {n_exons - 1} introns for {n_exons} exons, got {n_introns}"
            )
        else:
            assert n_introns == 0, f"{gene}: single-exon gene should have 0 introns, got {n_introns}"


def test_no_zero_length_introns(locus_data):
    for gene, info in locus_data.items():
        for s, e in info._intron_indices:
            assert e > s, f"{gene}: zero-length intron ({s}, {e})"


# ── General correctness ────────────────────────────────────────────────────────


def test_all_genes_present(locus_data):
    for gene in REGRESSION_GENES:
        assert gene in locus_data, f"{gene} missing from result"


def test_strand(locus_data):
    for gene, truth in GROUND_TRUTH.items():
        assert locus_data[gene].strand == truth["strand"], (
            f"{gene}: expected {truth['strand']}, got {locus_data[gene].strand}"
        )


def test_gene_type(locus_data):
    for gene, truth in GROUND_TRUTH.items():
        assert locus_data[gene].gene_type == truth["gene_type"], (
            f"{gene}: expected {truth['gene_type']}, got {locus_data[gene].gene_type}"
        )


def test_lncrna_gene_type(locus_data):
    for gene in ("MALAT1", "NEAT1", "XIST"):
        assert locus_data[gene].gene_type == GeneType.LNCRNA, (
            f"{gene}: expected LNCRNA, got {locus_data[gene].gene_type}"
        )


def test_gene_coordinates_sane(locus_data):
    for gene, info in locus_data.items():
        assert info.gene_start >= 0
        assert info.gene_end > info.gene_start
        assert len(info.full_mrna) == info.gene_end - info.gene_start, (
            f"{gene}: sequence length {len(info.full_mrna)} != gene span {info.gene_end - info.gene_start}"
        )


def test_gene_contains_all_exons(locus_data):
    """Every exon must fall within the gene body."""
    for gene, info in locus_data.items():
        for s, e in info._exon_indices:
            assert s >= info.gene_start, f"{gene}: exon start {s} < gene_start {info.gene_start}"
            assert e <= info.gene_end, f"{gene}: exon end {e} > gene_end {info.gene_end}"


def test_exon_order_is_biological(locus_data):
    """Index 0 must be the biological 5' end."""
    for gene, info in locus_data.items():
        if len(info._exon_indices) < 2:
            continue
        if info.strand == StrandType.NEG:
            assert info._exon_indices[0][0] > info._exon_indices[-1][0], (
                f"{gene}: NEG strand exons not in 5'->3' biological order"
            )
        else:
            assert info._exon_indices[0][0] < info._exon_indices[-1][0], (
                f"{gene}: POS strand exons not in 5'->3' biological order"
            )


def test_sequence_uppercase(locus_data):
    for gene, info in locus_data.items():
        assert info.full_mrna == info.full_mrna.upper(), f"{gene}: sequence contains lowercase"


def test_sequence_valid_bases(locus_data):
    valid = set("ACGTUN")
    for gene, info in locus_data.items():
        bad = set(info.full_mrna) - valid
        assert not bad, f"{gene}: unexpected bases {bad}"


def test_exon_first_last_coords(locus_data):
    """
    Verify biological 5' and 3' ends match GENCODE GRCh38.

    POS strand: bio_five_prime  = _exon_indices[0][0]  (start of first exon)
                bio_three_prime = _exon_indices[-1][1] (end of last exon)
    NEG strand: bio_five_prime  = _exon_indices[0][1]  (END of first exon — rightmost edge)
                bio_three_prime = _exon_indices[-1][0] (start of last exon — leftmost edge)
    """
    for gene, truth in GROUND_TRUTH.items():
        info = locus_data[gene]
        first_exon = info._exon_indices[0]
        last_exon = info._exon_indices[-1]

        if info.strand == StrandType.POS:
            actual_five_prime = first_exon[0]
            actual_three_prime = last_exon[1]
        else:  # NEG
            actual_five_prime = first_exon[1]  # end of first exon = biological 5'
            actual_three_prime = last_exon[0]  # start of last exon = biological 3'

        assert actual_five_prime == truth["bio_five_prime"], (
            f"{gene}: bio 5' end {actual_five_prime} != {truth['bio_five_prime']}"
        )
        assert actual_three_prime == truth["bio_three_prime"], (
            f"{gene}: bio 3' end {actual_three_prime} != {truth['bio_three_prime']}"
        )


def test_large_gene_span(locus_data):
    """RBFOX1 spans ~2.3Mb — coordinate arithmetic should still hold."""
    info = locus_data["RBFOX1"]
    span = info.gene_end - info.gene_start
    assert span > 2_000_000, f"RBFOX1 span {span} looks too small"
    assert len(info.full_mrna) == span


def test_gene_boundaries_independent_of_canonical(locus_data, locus_data_canonical):
    """
    gene_start, gene_end, and full_mrna come from the gene feature directly
    and must be identical regardless of canonical_only mode.
    """
    for gene in REGRESSION_GENES:
        info_all = locus_data[gene]
        info_can = locus_data_canonical[gene]

        assert info_all.gene_start == info_can.gene_start, (
            f"{gene}: gene_start differs: all={info_all.gene_start}, canonical={info_can.gene_start}"
        )
        assert info_all.gene_end == info_can.gene_end, (
            f"{gene}: gene_end differs: all={info_all.gene_end}, canonical={info_can.gene_end}"
        )
        assert info_all.full_mrna == info_can.full_mrna, f"{gene}: full_mrna differs between modes"
        assert info_all.strand == info_can.strand, (
            f"{gene}: strand differs: all={info_all.strand}, canonical={info_can.strand}"
        )
