"""Cross-feature interaction features."""

from ..common.modifications import deoxy_sugar_fraction, phosphorothioate_fraction

# Expected strongest self-folding, about -0.5 kcal/mol per nucleotide; normalizing the self-fold MFE
# by it puts the structure term on a unit scale comparable to the [0,1] sugar/backbone fractions.
STACK_REF_KCAL_PER_NT = 0.5


def protein_affinity_gymnosis(ps_pattern, chemical_pattern, seq_internal_fold, transfection_gymnosis):
    """ASO -> cell-surface-protein affinity, gated to gymnotic uptake (Gaus et al. 2019, NAR 47(3):1110).

    Sum of the PS fraction, the deoxy-sugar fraction and the self-folding MFE -- the last normalized
    according to the expected strongest folding of about -0.5 kcal/mol per nucleotide -- multiplied by
    the gymnosis factor. All inputs are aligned pandas Series (gymnosis is 0 for assisted/unknown delivery).
    """
    length = chemical_pattern.astype(str).str.len().clip(lower=1)
    ps_frac = ps_pattern.map(phosphorothioate_fraction)
    dna_frac = chemical_pattern.map(deoxy_sugar_fraction)
    affinity = ps_frac + dna_frac + seq_internal_fold / (STACK_REF_KCAL_PER_NT * length)
    is_gymnosis = transfection_gymnosis.fillna(0.0) if transfection_gymnosis is not None else 0.0
    return affinity * is_gymnosis
