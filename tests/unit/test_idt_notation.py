"""IDT order-string rendering for ASO chemistries (PS-DNA, 2'-MOE, cEt, LNA)."""

import pytest

from tauso.common.modifications import to_idt_notation


def test_ps_dna_full_backbone():
    # plain DNA, all phosphorothioate -> bases joined by '*'
    assert to_idt_notation("ACGT", "dddd", "***") == "A*C*G*T"


def test_ps_pattern_defaults_to_all_ps():
    assert to_idt_notation("ACGT", "dddd") == "A*C*G*T"


def test_phosphodiester_linkage_emits_no_star():
    # middle linkage phosphodiester ('o') -> no '*' between C and G
    assert to_idt_notation("ACGT", "dddd", "*o*") == "A*CG*T"


def test_moe_gapmer_terminal_and_internal_codes():
    # 2-2-2 MOE gapmer: 5'/3' termini use terminal codes, inner wing residue uses the internal code
    out = to_idt_notation("AACGTT", "MMddMM", "*****")
    assert out == "/52MOErA/*/i2MOErA/*C*G*/i2MOErT/*/32MOErT/"


def test_lna_gapmer_uses_plus_prefix():
    out = to_idt_notation("AACGTT", "LLddLL", "*****")
    assert out == "+A*+A*C*G*+T*+T"


def test_cet_gapmer_uses_cet_code():
    # pattern letter 'C' is the cEt sugar; the DNA-gap base 'C' stays plain
    out = to_idt_notation("AACGTT", "CCddCC", "*****")
    assert out == "/icEtA/*/icEtA/*C*G*/icEtT/*/icEtT/"


@pytest.mark.parametrize(
    "seq, pattern, ps",
    [
        ("ACGT", "ddd", "***"),   # chemical_pattern shorter than sequence
        ("ACGT", "dddd", "**"),   # ps_pattern is not len(seq)-1
        ("ACGT", "ddXd", "***"),  # unsupported sugar 'X'
        ("ACGZ", "dddd", "***"),  # unsupported base 'Z'
    ],
)
def test_invalid_inputs_raise(seq, pattern, ps):
    with pytest.raises(ValueError):
        to_idt_notation(seq, pattern, ps)
