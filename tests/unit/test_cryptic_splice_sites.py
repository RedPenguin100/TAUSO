from types import SimpleNamespace

import numpy as np
import pytest

from tauso.genome.LocusInfo import LocusInfo, StrandType
from tauso.populate.populate_structure import (
    _SPLICE_SCAN_WINDOW,
    assign_cryptic_splice_sites,
    maxent_scorers,
)

DONOR = "CAGGTAAGT"
GENE_LENGTH = 600
INTRON = (200, 400)


def locus_with(sequence, strand=StrandType.POS):
    locus = LocusInfo()
    locus.full_mrna = sequence
    locus.gene_start, locus.gene_end, locus.strand = 0, len(sequence), strand
    locus.add_exon_indices(0, INTRON[0])
    locus.add_intron_indices(*INTRON)
    locus.add_exon_indices(INTRON[1], len(sequence))
    return locus


def gene_sequence(cryptic_at=None):
    """A gene whose intron carries a real donor, optionally with a second donor planted inside."""
    seq = list("GCAT" * (GENE_LENGTH // 4))
    seq[INTRON[0] - 3 : INTRON[0] + 6] = list(DONOR)
    if cryptic_at is not None:
        seq[cryptic_at : cryptic_at + len(DONOR)] = list(DONOR)
    return "".join(seq)


def run(locus, coordinates):
    out = SimpleNamespace(
        **{
            name: np.full(len(coordinates), np.nan)
            for name in (
                "cryptic_donor_max",
                "cryptic_acceptor_max",
                "local_donor_mean",
                "local_acceptor_mean",
                "cryptic_donor_delta",
                "cryptic_acceptor_delta",
            )
        }
    )
    assign_cryptic_splice_sites(out, np.arange(len(coordinates)), np.array(coordinates), locus, maxent_scorers())
    return out


def test_planted_donor_is_found_in_the_window():
    target = 300
    out = run(locus_with(gene_sequence(cryptic_at=target)), [target])
    assert out.cryptic_donor_max[0] == pytest.approx(10.858, abs=0.01)


def test_window_scores_are_defined_without_a_host_intron():
    out = run(locus_with(gene_sequence()), [100])  # exonic
    assert np.isfinite(out.cryptic_donor_max[0])
    assert np.isfinite(out.local_donor_mean[0])
    assert np.isnan(out.cryptic_donor_delta[0])


def test_delta_is_defined_only_inside_an_intron():
    sequence = gene_sequence()
    out = run(locus_with(sequence), [100, 300])
    assert np.isnan(out.cryptic_donor_delta[0])
    assert np.isfinite(out.cryptic_donor_delta[1])


def test_a_cryptic_donor_matching_the_real_one_gives_delta_zero():
    target = 300
    out = run(locus_with(gene_sequence(cryptic_at=target)), [target])
    assert out.cryptic_donor_delta[0] == pytest.approx(0.0, abs=0.01)


def test_plain_intron_scores_below_its_own_donor():
    out = run(locus_with(gene_sequence()), [300])
    assert out.cryptic_donor_delta[0] < 0


def test_max_is_at_least_the_mean():
    out = run(locus_with(gene_sequence(cryptic_at=300)), [300])
    assert out.cryptic_donor_max[0] >= out.local_donor_mean[0]
    assert out.cryptic_acceptor_max[0] >= out.local_acceptor_mean[0]


def test_donor_outside_the_window_is_not_picked_up():
    target = 300
    far = target + _SPLICE_SCAN_WINDOW + 20
    near = run(locus_with(gene_sequence(cryptic_at=target)), [target])
    away = run(locus_with(gene_sequence(cryptic_at=far)), [target])
    assert away.cryptic_donor_max[0] < near.cryptic_donor_max[0]


def test_rna_and_dna_score_alike():
    sequence = gene_sequence(cryptic_at=300)
    dna = run(locus_with(sequence), [300])
    rna = run(locus_with(sequence.replace("T", "U")), [300])
    assert dna.cryptic_donor_max[0] == rna.cryptic_donor_max[0]
    assert dna.cryptic_donor_delta[0] == rna.cryptic_donor_delta[0]


def test_ambiguous_bases_do_not_raise():
    sequence = list(gene_sequence(cryptic_at=300))
    sequence[310:330] = list("N" * 20)
    out = run(locus_with("".join(sequence)), [300])
    assert np.isfinite(out.cryptic_donor_max[0])


def test_an_all_n_window_is_nan():
    sequence = list(gene_sequence())
    sequence[240:360] = list("N" * 120)
    out = run(locus_with("".join(sequence)), [300])
    assert np.isnan(out.cryptic_donor_max[0])
    assert np.isnan(out.local_donor_mean[0])
