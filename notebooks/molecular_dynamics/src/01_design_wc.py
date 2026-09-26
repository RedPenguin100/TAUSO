#!/usr/bin/env python
"""Step 1: the uniform arm -- one sequence set, four chemistries, each against RNA.

Every line of `config/design/design.json` is one duplex: the two strand sequences, and which sugar
sits at each position. Everything downstream builds and runs what is on this list.

Only duplexes against an RNA strand are designed here -- the state an antisense oligo is in when
bound to its target. A 5-10-5 gapmer on target is `M:R` in its wings, `D:R` in its gap, and
`MD|RR` / `DM|RR` at the two boundaries between them.

    48 sequences x 4 sugars against ribose = 192 duplexes

The sequences are generated here, under `usable`, and written to `config/design/sequences.json`.
Every arm then reads that one file, so a cell measured in one chemistry and a cell measured in
another describe the same sequences. Under the sixteen-cell split each uniform cell gets 32-34
steps drawn from 25-32 duplexes, so no extra sequence design is needed for this arm; the junction
arm draws its own in step 2.

Coverage is counted in duplexes, not steps. Replicates of one duplex differ only in starting
velocities, so a cell measured in two duplexes cannot be given an uncertainty by resampling
however many replicates it has -- three duplexes is the minimum for an error bar.
"""
import json

import numpy as np

from consts import DESIGN, DUPLEX_LENGTH as L, describe_design, ensure_dirs, reverse_complement
from sequence_design import (CANDIDATE_POOL, N_SEQUENCES, SEQUENCE_SEED, balanced, make_sequences)

SEQUENCES = DESIGN.parent / "sequences.json"                 # written here, then read by every arm

# the sugar on strand A; strand B is ribose throughout, so each state reads "<sugar>:R"
#   D  2'-deoxy, the gap of a gapmer        M  2'-MOE, a wing
#   E  cEt, a wing                          R  ribose, the RNA:RNA reference
SUGARS = "DMER"


def sequences():
    """The Watson-Crick sequences every arm is built from, generated here and then written down.

    Every sequence is drawn under `usable` and none is modified afterwards, so the file this
    writes is reproducible from the seed alone. Delete it to redraw.
    """
    if SEQUENCES.exists():
        return json.loads(SEQUENCES.read_text())["sequences"]
    pool = make_sequences(CANDIDATE_POOL, np.random.default_rng(SEQUENCE_SEED))
    seqs = balanced(pool, N_SEQUENCES)
    SEQUENCES.write_text(json.dumps({
        "note": (f"{N_SEQUENCES} Watson-Crick sequences of length {L}, drawn under usable() with "
                 f"numpy default_rng({SEQUENCE_SEED}): a pool of {CANDIDATE_POOL} drawn under "
                 f"usable(), thinned to the subset covering the sixteen dinucleotides most "
                 f"evenly. Every one passes the filter and none was modified after selection, so "
                 f"this file regenerates from the seed alone. Delete it to redraw."),
        "seed": SEQUENCE_SEED, "length": L, "sequences": seqs}, indent=1))
    return seqs


def main():
    ensure_dirs()
    rows = [dict(name=f"R{sugar}_{i:03d}", family="rna", state=f"{sugar}:R",
                 status="todo", seq=seq, partner=reverse_complement(seq),
                 sugarA=sugar, sugarB="R")
            for sugar in SUGARS
            for i, seq in enumerate(sequences())]

    DESIGN.write_text(json.dumps(rows, indent=1))
    describe_design(rows)


if __name__ == "__main__":
    main()
