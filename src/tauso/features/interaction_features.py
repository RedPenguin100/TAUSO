"""Cross-feature interaction features."""


def internal_fold_gymnosis(self_fold, transfection_gymnosis):
    """ASO self-fold gated to gymnotic (carrier-free) uptake.

    A self-structured ASO is less available to be taken up and to hybridise when the cell must take it
    up unaided (gymnosis); a transfection reagent or electroporation bypasses this. The value is the
    self-fold energy ``self_fold`` (self-fold MFE, <= 0) on gymnosis rows and 0
    elsewhere. Inputs are aligned pandas Series; ``transfection_gymnosis`` is the 0/1 gymnosis gate
    (None -> the feature is 0 for all rows).
    """
    if transfection_gymnosis is None:
        return self_fold * 0.0
    is_gymnosis = transfection_gymnosis.fillna(0.0)
    return self_fold.where(is_gymnosis > 0, 0.0)
