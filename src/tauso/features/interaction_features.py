"""Cross-feature interaction features."""


def internal_fold_gymnosis(internal_fold, transfection_gymnosis):
    """ASO internal-fold energy (MFE, <= 0) gated to gymnotic (carrier-free) uptake: the fold energy on
    gymnosis rows, 0 elsewhere. A self-structured ASO is less available for unaided (gymnotic) uptake and
    hybridisation; a transfection reagent or electroporation bypasses this. Inputs are aligned Series."""
    return internal_fold.where(transfection_gymnosis.fillna(0.0) > 0, 0.0)
