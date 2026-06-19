"""Cross-feature interaction features."""

from ..common.modifications import deoxy_sugar_fraction, phosphorothioate_fraction

# Scaling constant (kcal/mol per nucleotide) that renders the per-nucleotide self-folding energy
# dimensionless and comparable to the [0,1] sugar/backbone fractions. Scale-setting, not fitted: the
# feature is beneficial across kappa in [0.125, 0.5] (a 4x range), so its exact value is not critical.
STACK_REF_KCAL_PER_NT = 0.25


def protein_affinity_gymnosis(ps_pattern, chemical_pattern, seq_internal_fold, transfection_gymnosis):
    """ASO -> cell-surface-protein affinity, gated to gymnotic uptake (Gaus et al. 2019, NAR 47(3):1110).

    Sum of the PS fraction, the deoxy-sugar fraction and the self-folding MFE (seq_internal_fold,
    the DNA-parameter self-fold), the last normalized per nucleotide by STACK_REF_KCAL_PER_NT, all
    multiplied by the gymnosis factor. Inputs are aligned pandas Series (gymnosis 0 for assisted/unknown).
    """
    length = chemical_pattern.astype(str).str.len().clip(lower=1)
    ps_frac = ps_pattern.map(phosphorothioate_fraction)
    dna_frac = chemical_pattern.map(deoxy_sugar_fraction)
    affinity = ps_frac + dna_frac + seq_internal_fold / (STACK_REF_KCAL_PER_NT * length)
    is_gymnosis = transfection_gymnosis.fillna(0.0) if transfection_gymnosis is not None else 0.0
    return affinity * is_gymnosis
