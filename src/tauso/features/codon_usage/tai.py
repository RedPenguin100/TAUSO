import os
from dataclasses import dataclass
from enum import Enum

import numpy as np
import pandas as pd
from codonbias.scores import TrnaAdaptationIndex

from ...data.data import get_data_dir


@dataclass(frozen=True)
class TGCNSourceSpec:
    """Where to fetch a tGCN table from and how to verify it. One entry per
    specific GtRNAdb genome (e.g. Hsapi38 for human, Mmusc10 for mouse)."""

    organism: str  # display name, e.g. "human"
    filename: str  # filename inside TAUSO_DATA_DIR
    gtrnadb_genome: str  # GtRNAdb genome ID, e.g. "Hsapi38"
    gtrnadb_domain: str  # GtRNAdb taxonomic domain, e.g. "eukaryota"
    sha256: str  # pinned SHA256 of the on-disk CSV


class TGCNSource(Enum):
    """Supported tGCN sources, keyed by specific GtRNAdb genome ID."""

    HSAPI38 = TGCNSourceSpec(
        organism="human",
        filename="human_tgcn_hsapi38.csv",
        gtrnadb_genome="Hsapi38",
        gtrnadb_domain="eukaryota",
        sha256="80bbf6c62395e6b9463e56db16efa87fccdc3c6d2452b808099482c33f23d8e6",
    )
    # MMUSC10 = TGCNSourceSpec(
    #     organism="mouse",
    #     filename="mouse_tgcn_mmusc10.csv",
    #     gtrnadb_genome="Mmusc10",
    #     gtrnadb_domain="eukaryota",
    #     sha256="<pin>",
    # )


_SCORER_CACHE: dict[TGCNSource, TrnaAdaptationIndex] = {}


def _load_scorer(source: TGCNSource) -> TrnaAdaptationIndex:
    cached = _SCORER_CACHE.get(source)
    if cached is not None:
        return cached
    spec = source.value
    path = os.path.join(get_data_dir(), spec.filename)
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Missing tGCN file at {path}. Run 'tauso setup-tgcn --source {source.name.lower()}' to download it."
        )
    tGCN = pd.read_csv(path)
    scorer = TrnaAdaptationIndex(tGCN=tGCN, s_values="dosReis", genetic_code=1)
    _SCORER_CACHE[source] = scorer
    return scorer


def compute_tAI(seq: str, source: TGCNSource = TGCNSource.HSAPI38) -> float:
    if not isinstance(seq, str) or not seq.strip():
        return np.nan
    score = _load_scorer(source).get_score(seq)
    return float(score) if np.isfinite(score) else np.nan
