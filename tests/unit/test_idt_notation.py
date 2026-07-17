"""IDT order-string rendering for ASO chemistries (PS-DNA, 2'-MOE, LNA, 5-methyl-C)."""

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


def test_cet_sugar_is_rejected():
    # cEt (pattern letter 'C') is not an IDT catalogue modification -> ValueError, not a rendered code
    with pytest.raises(ValueError, match="cEt"):
        to_idt_notation("AACGTT", "CCddCC", "*****")


def test_moe_gapmer_5methyl_c_gap():
    # 5-methylcytosine chemistry: DNA-gap C -> /iMe-dC/; gap A/G/T unchanged
    out = to_idt_notation("AACGTT", "MMddMM", "*****", "MOE/5-methylcytosines/deoxy")
    assert out == "/52MOErA/*/i2MOErA/*/iMe-dC/*G*/i2MOErT/*/32MOErT/"


def test_5methyl_c_wing_c_already_methyl_via_moe_code():
    # a wing C stays /..2MOErC/ (2'-MOE C is already 5-methyl); only the DNA-gap C becomes /iMe-dC/
    out = to_idt_notation("CACGTC", "MMddMM", "*****", "MOE/5-methylcytosines/deoxy")
    assert out == "/52MOErC/*/i2MOErA/*/iMe-dC/*G*/i2MOErT/*/32MOErC/"


def test_no_modification_string_leaves_gap_c_plain():
    # without the 5-methylcytosine label, DNA-gap C stays plain 'C' (backwards compatible)
    assert to_idt_notation("AACGTT", "MMddMM", "*****") == "/52MOErA/*/i2MOErA/*C*G*/i2MOErT/*/32MOErT/"


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
