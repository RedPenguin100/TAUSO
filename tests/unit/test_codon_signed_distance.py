"""Sign convention of the signed start/stop-codon distances.

The gene is synthetic so the nearest codon of every coordinate is known by hand: a start
codon at [100, 103) and a stop codon at [200, 203), whose midpoints are 101.5 and 201.5.
"""
import numpy as np
import pytest

from tauso.genome.LocusInfo import LocusInfo, StrandType
from tauso.populate.populate_structure import _init_outputs, assign_distance_to_special_codons

START, STOP = (100, 103), (200, 203)
START_MID, STOP_MID = 101.5, 201.5


def locus(strand, all_starts=None, all_stops=None):
    info = LocusInfo()
    info.start_codon, info.stop_codon = START, STOP
    info.all_start_codons = [START] if all_starts is None else all_starts
    info.all_stop_codons = [STOP] if all_stops is None else all_stops
    info.gene_start, info.gene_end = 0, 300
    info.strand = strand
    return info


def run(coords, strand, **kwargs):
    coords = np.asarray(coords, dtype=np.float64)
    out = _init_outputs(len(coords))
    assign_distance_to_special_codons(out, np.arange(len(coords)), coords, locus(strand, **kwargs))
    return out


@pytest.mark.parametrize("coord, codon_is_5_prime", [(150, True), (50, False)])
def test_sign_follows_codon_direction_on_plus_strand(coord, codon_is_5_prime):
    """Positive when the codon lies 5' of the target's centre."""
    out = run([coord], StrandType.POS)
    expected = abs(coord - START_MID) * (1 if codon_is_5_prime else -1)
    assert out.signed_dist_canonical_start[0] == pytest.approx(expected)


@pytest.mark.parametrize("coord, codon_is_5_prime", [(150, False), (50, True)])
def test_minus_strand_flips_the_sign(coord, codon_is_5_prime):
    """Transcript orientation, not genomic: the same coordinate reverses on the minus strand."""
    out = run([coord], StrandType.NEG)
    expected = abs(coord - START_MID) * (1 if codon_is_5_prime else -1)
    assert out.signed_dist_canonical_start[0] == pytest.approx(expected)


def test_magnitude_matches_the_unsigned_distance():
    coords = [0, 50, 101, 150, 202, 299]
    for strand in (StrandType.POS, StrandType.NEG):
        out = run(coords, strand)
        for signed, unsigned in (
            (out.signed_dist_canonical_start, out.dist_canonical_start),
            (out.signed_dist_canonical_stop, out.dist_canonical_stop),
            (out.signed_dist_closest_start, out.dist_closest_start),
            (out.signed_dist_closest_stop, out.dist_closest_stop),
        ):
            np.testing.assert_allclose(np.abs(signed), unsigned)


def test_start_and_stop_are_measured_against_their_own_codon():
    """A coordinate between the two is 5' of the stop and 3' of the start."""
    out = run([150], StrandType.POS)
    assert out.signed_dist_canonical_start[0] == pytest.approx(150 - START_MID)  # start is 5', positive
    assert out.signed_dist_canonical_stop[0] == pytest.approx(150 - STOP_MID)  # stop is 3', negative
    assert out.signed_dist_canonical_stop[0] < 0


def test_closest_picks_the_nearest_of_several_isoform_codons():
    out = run([150], StrandType.POS, all_starts=[(100, 103), (140, 143)])
    assert out.signed_dist_closest_start[0] == pytest.approx(150 - 141.5)
    assert out.signed_dist_canonical_start[0] == pytest.approx(150 - START_MID)


def test_non_coding_gene_stays_nan():
    info = LocusInfo()
    info.start_codon = info.stop_codon = None
    info.all_start_codons = info.all_stop_codons = []
    info.gene_start, info.gene_end = 0, 300
    info.strand = StrandType.POS
    out = _init_outputs(1)
    assign_distance_to_special_codons(out, np.arange(1), np.array([150.0]), info)
    for arr in (out.signed_dist_canonical_start, out.signed_dist_canonical_stop,
                out.signed_dist_closest_start, out.signed_dist_closest_stop):
        assert np.isnan(arr[0])
