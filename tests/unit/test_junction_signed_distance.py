"""Sign convention and side-filing of the canonical splice-junction distances.

The gene is synthetic so the nearest junction of every coordinate is known by hand:
exons [0, 100) and [200, 300) with the single intron [100, 200) between them, which
puts junctions at 100 and 200.
"""
import numpy as np
import pytest

from tauso.genome.LocusInfo import LocusInfo, StrandType
from tauso.populate.populate_structure import (
    _init_outputs,
    assign_canonical_splice_junction_distances,
)

EXONS = [(0, 100), (200, 300)]
INTRONS = [(100, 200)]


def locus(strand):
    info = LocusInfo()
    info._exon_indices = EXONS
    info._intron_indices = INTRONS
    info.gene_start, info.gene_end = 0, 300
    info.strand = strand
    return info


def run(coords, strand):
    """Assign the junction features for `coords` and return the outputs object."""
    coords = np.asarray(coords)
    out = _init_outputs(len(coords))
    assign_canonical_splice_junction_distances(out, np.arange(len(coords)), coords, locus(strand))
    return out


@pytest.mark.parametrize(
    "coord, distance, junction_is_5_prime",
    [
        (50, 50, False),  # in the first exon, junction 100 lies 3'
        (90, 10, False),
        (250, 50, True),  # in the second exon, junction 200 lies 5'
        (210, 10, True),
    ],
)
def test_sign_follows_junction_direction_on_plus_strand(coord, distance, junction_is_5_prime):
    out = run([coord], StrandType.POS)
    expected = np.log1p(distance) * (1 if junction_is_5_prime else -1)
    assert out.junction_signed_logdist_exonic[0] == pytest.approx(expected)
    assert np.isnan(out.junction_signed_logdist_intronic[0])


@pytest.mark.parametrize("coord, distance, junction_is_5_prime", [(50, 50, True), (250, 50, False)])
def test_minus_strand_flips_the_sign(coord, distance, junction_is_5_prime):
    """Transcript orientation, not genomic: the same coordinate reverses on the minus strand."""
    out = run([coord], StrandType.NEG)
    expected = np.log1p(distance) * (1 if junction_is_5_prime else -1)
    assert out.junction_signed_logdist_exonic[0] == pytest.approx(expected)


def test_intronic_site_is_filed_on_the_intronic_side_only():
    out = run([150], StrandType.POS)
    assert np.isnan(out.junction_signed_logdist_exonic[0])
    assert out.junction_signed_logdist_intronic[0] == pytest.approx(np.log1p(50))


def test_site_in_neither_exon_nor_intron_stays_nan_on_both_sides():
    """A coordinate past the last exon is in neither, so no side may claim it."""
    out = run([500], StrandType.POS)
    assert np.isnan(out.junction_signed_logdist_exonic[0])
    assert np.isnan(out.junction_signed_logdist_intronic[0])
    assert np.isnan(out.dist_sj_exonic[0])
    assert np.isnan(out.dist_sj_intronic[0])


def test_magnitude_matches_the_unsigned_log_distance():
    coords = [50, 90, 150, 210, 250]
    out = run(coords, StrandType.POS)
    signed = np.where(np.isnan(out.junction_signed_logdist_exonic),
                      out.junction_signed_logdist_intronic, out.junction_signed_logdist_exonic)
    unsigned = np.where(np.isnan(out.junction_logdist_exonic),
                        out.junction_logdist_intronic, out.junction_logdist_exonic)
    np.testing.assert_allclose(np.abs(signed), unsigned)


def test_single_exon_gene_leaves_every_junction_feature_nan():
    info = LocusInfo()
    info._exon_indices = [(0, 300)]
    info._intron_indices = []
    info.gene_start, info.gene_end = 0, 300
    info.strand = StrandType.POS
    out = _init_outputs(1)
    assign_canonical_splice_junction_distances(out, np.arange(1), np.array([50]), info)
    assert np.isnan(out.junction_signed_logdist_exonic[0])
    assert np.isnan(out.junction_signed_logdist_intronic[0])
