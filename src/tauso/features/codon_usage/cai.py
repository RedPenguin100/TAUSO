from __future__ import annotations

import numpy as np
import pandas as pd
from codonbias.scores import CodonAdaptationIndex

# JSON written by `tauso build-cai-weights`: {cell_line: {codon: weight}}
# plus a "Generic" entry holding the cohort-consensus fallback.
CAI_WEIGHTS_FILENAME = "cai_weights.json"

# Default pseudocount for the build-time weight computation. 0 matches the
# previous hand-rolled `calc_CAI_weight` (raw codon counts normalised per
# amino acid by the max, no smoothing). codonbias's own default is 1, which
# avoids zero weights for unobserved codons at the cost of a small shift in
# every weight; callers can override via the CLI.
CAI_DEFAULT_PSEUDOCOUNT = 0


def _all_codons_seq() -> str:
    """Return a synthetic CDS containing every one of the 64 codons exactly
    once. Used to bootstrap a CodonAdaptationIndex with non-zero counts for
    every codon, so codonbias's internal log() never sees a zero (avoids
    spurious RuntimeWarnings) before we overwrite the weights anyway.
    """
    bases = "ACGT"
    return "".join(b1 + b2 + b3 for b1 in bases for b2 in bases for b3 in bases)


_SCORER_BOOTSTRAP_SEQ = _all_codons_seq()


def _make_scorer_from_weights(weights: dict[str, float]) -> CodonAdaptationIndex:
    """Construct a CodonAdaptationIndex and overwrite its internal weight
    tables with `weights`. The bootstrap sequence above is only there to
    let __init__ build `counter`, `kmer_index`, etc.; we then replace the
    learned weights with the saved ones.

    TODO(codonbias-internals): this writes to attributes that are not part
    of codonbias's public API (`weights`, `log_weights`, `_log_weights_arr`).
    Replace with a public hook once codon-bias upstream exposes one.
    """
    scorer = CodonAdaptationIndex(ref_seq=[_SCORER_BOOTSTRAP_SEQ], pseudocount=CAI_DEFAULT_PSEUDOCOUNT)
    scorer.weights = pd.Series(weights, name="weight")
    scorer.log_weights = np.log(scorer.weights)
    scorer._log_weights_arr = scorer.log_weights.reindex(scorer.counter.kmer_index).values
    return scorer


def build_scorer_from_reference(
    cds_list: list[str], pseudocount: int = CAI_DEFAULT_PSEUDOCOUNT
) -> CodonAdaptationIndex:
    """Build a CodonAdaptationIndex from a list of reference CDS sequences.
    Used at build-time by `tauso build-cai-weights`. See codonbias's own
    docs for the meaning of `pseudocount`."""
    return CodonAdaptationIndex(ref_seq=cds_list, pseudocount=pseudocount)


class CAIScorerCache:
    """Lazy per-cell-line CAI scorer cache. Wraps codonbias scorers whose
    weight tables are loaded from the saved cai_weights.json (one weight
    dict per cell line, plus a 'Generic' fallback)."""

    GENERIC_KEY = "Generic"

    def __init__(self, weights_map: dict[str, dict[str, float]]):
        self._weights_map = weights_map
        self._scorers: dict[str, CodonAdaptationIndex | None] = {}

    def _build_scorer(self, cell_line: str) -> CodonAdaptationIndex | None:
        weights = self._weights_map.get(cell_line) or self._weights_map.get(self.GENERIC_KEY)
        if not weights:
            return None
        return _make_scorer_from_weights(weights)

    def get_scorer(self, cell_line: str) -> CodonAdaptationIndex | None:
        if cell_line not in self._scorers:
            self._scorers[cell_line] = self._build_scorer(cell_line)
        return self._scorers[cell_line]

    def score(self, seq: str, cell_line: str) -> float:
        if not isinstance(seq, str) or not seq.strip():
            return np.nan
        scorer = self.get_scorer(cell_line)
        if scorer is None:
            return np.nan
        value = scorer.get_score(seq)
        return float(value) if np.isfinite(value) else np.nan

    def known_cell_lines(self) -> set[str]:
        return set(self._weights_map.keys())
