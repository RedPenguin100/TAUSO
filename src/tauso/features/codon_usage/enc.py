import numpy as np
from codonbias.scores import EffectiveNumberOfCodons

# Match the parameters of our previous custom numba implementation:
# pseudocount=1, robust chi2 formula, background correction from BNC, weighted mean across AAs.
_ENC_SCORER = EffectiveNumberOfCodons(
    k_mer=1,
    bg_correction=True,
    robust=True,
    pseudocount=1,
    mean="weighted",
    genetic_code=1,
)


def compute_ENC(seq: str) -> float:
    if not isinstance(seq, str) or not seq.strip():
        return np.nan

    enc_score = _ENC_SCORER.get_score(seq)

    if np.isnan(enc_score):
        return np.nan

    return (enc_score - 20.0) / (61.0 - 20.0)
