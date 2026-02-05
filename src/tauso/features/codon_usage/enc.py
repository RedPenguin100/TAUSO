import numpy as np
import codonbias as cb

# TODO: this is way too slow for what its doing
def compute_ENC(seq: str) -> float:
    """
    Returns normalized ENC in [0,1], or NaN if the input sequence is empty/invalid.
    """
    if not isinstance(seq, str) or seq.strip() == "":
        return np.nan

    enc = cb.scores.EffectiveNumberOfCodons(bg_correction=True)
    enc_score = enc.get_score(seq)

    if enc_score is None or not np.isfinite(enc_score):
        return np.nan

    # Normalize ENC from [20..61] to [0..1]
    return (enc_score - 20.0) / (61.0 - 20.0)
