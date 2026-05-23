import os

import numpy as np
import pandas as pd
from codonbias.scores import TrnaAdaptationIndex

from ...data.data import get_data_dir

# Filename of the human tGCN snapshot written by `tauso setup-tgcn`.
# Source: GtRNAdb Hsapi38 (eukaryota), sorted by anti_codon.
HUMAN_TGCN_FILENAME = "human_tgcn_hsapi38.csv"

_TAI_SCORER = None


def _load_scorer() -> TrnaAdaptationIndex:
    path = os.path.join(get_data_dir(), HUMAN_TGCN_FILENAME)
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Missing tGCN file at {path}. Run 'tauso setup-tgcn' (or 'tauso setup-omics') to download it from GtRNAdb."
        )
    tGCN = pd.read_csv(path)
    return TrnaAdaptationIndex(tGCN=tGCN, s_values="dosReis", genetic_code=1)


def compute_tAI(seq: str) -> float:
    global _TAI_SCORER
    if _TAI_SCORER is None:
        _TAI_SCORER = _load_scorer()
    if not isinstance(seq, str) or not seq.strip():
        return np.nan
    score = _TAI_SCORER.get_score(seq)
    return float(score) if np.isfinite(score) else np.nan
