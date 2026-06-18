"""Cross-feature interaction features."""

from ..common.modifications import deoxy_sugar_fraction, phosphorothioate_fraction

# component stds on the relaxed (no-mixmer) population -> comparable scale before summing
PS_FRAC_STD, DNA_FRAC_STD, FOLD_STD = 0.17672, 0.09011, 1.52807


def protein_affinity_gymnosis(ps_pattern, chemical_pattern, seq_internal_fold, transfection_gymnosis):
    """ASO -> cell-surface-protein affinity, gated to gymnotic uptake (Gaus et al. 2019, NAR 47(3):1110).

    Std-scaled sum of PS fraction + deoxy-sugar fraction + self-folding MFE, times the gymnosis
    indicator (0 for assisted/unknown delivery). Aligned pandas Series in, Series out.
    """
    ps_frac = ps_pattern.map(phosphorothioate_fraction)
    dna_frac = chemical_pattern.map(deoxy_sugar_fraction)
    affinity = ps_frac / PS_FRAC_STD + dna_frac / DNA_FRAC_STD + seq_internal_fold / FOLD_STD
    is_gymnosis = transfection_gymnosis.fillna(0.0) if transfection_gymnosis is not None else 0.0
    return affinity * is_gymnosis
