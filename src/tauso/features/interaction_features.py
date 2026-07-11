"""Cross-feature interaction features."""


def internal_fold_gymnosis(internal_fold, transfection_gymnosis):
    """ASO internal-fold energy (MFE, <= 0) on gymnosis rows, 0 elsewhere. A self-structured ASO is less
    available for gymnotic uptake; a transfection reagent or electroporation bypasses this."""
    return internal_fold.where(transfection_gymnosis.fillna(0.0) > 0, 0.0)
