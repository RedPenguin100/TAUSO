#!/usr/bin/env python
"""Step 2: the junction arm -- a gapmer strand whose sugar changes, paired with ribose.

Against RNA, strand B is uniformly ribose, so only strand A carries a junction and there are two
patterns -- `MD|RR` where the wing meets the gap, and `DM|RR` where the gap meets the other wing.
A gapmer contains exactly one of each, so one duplex measures both.

Layout is a 5-6-5: positions 1-5 wing, 6-11 gap, 12-16 wing, putting the two junction steps at 5
and 11, both clear of the ends by the margin the campaign fits on.

Six duplexes per cell. The junction table is the one nothing published supplies and the one whose
cells are thinnest, so this is where extra sampling is worth spending: it halves the error on a
junction weight for 96 more duplexes.

The junction dinucleotides are written into a sequence after it is drawn, so every sequence is put
back through `usable` afterwards and redrawn until it passes.
"""
import json
from collections import Counter

import numpy as np

from consts import DESIGN, describe_design, ensure_dirs, reverse_complement
from sequence_design import make_sequences, usable

JUNCTION_RNA = DESIGN.parent / "junction_rna.json"
SITES = (5, 11)                 # 1-based step index of each junction
DUPLEXES_PER_CELL = 6
RETRIES = 400                   # redraws allowed before a stamped sequence is given up on
DINUCS = [a + b for a in "AGTC" for b in "AGTC"]
SEEDS = {"M": 17, "E": 23}      # one sugar per arm, each with its own stream


def sugar_pattern(modified):
    """Strand A: wing / gap / wing, so step 5 reads <mod>D and step 11 reads D<mod>."""
    return modified * 5 + "D" * 6 + modified * 5


def with_junctions(sequence, first, second):
    """Write the two junction dinucleotides into a sequence at `SITES`."""
    seq = list(sequence)
    for site, dinuc in zip(SITES, (first, second)):
        seq[site - 1], seq[site] = dinuc
    return "".join(seq)


def designs(modified, rng):
    """One duplex per (5' junction dinucleotide, 3' junction dinucleotide) assignment."""
    sugars = sugar_pattern(modified)
    out = []
    for rep in range(DUPLEXES_PER_CELL):
        for k in range(16):
            # each repeat pairs the sixteen 5' cells against a different rotation of the 3' ones,
            # so both ends of the gap see all sixteen without repeating a combination
            first, second = DINUCS[k], DINUCS[(k + 6 * rep) % 16]
            for _ in range(RETRIES):
                seq = with_junctions(make_sequences(1, rng)[0], first, second)
                if usable(seq):
                    break
            else:
                raise SystemExit(f"no usable sequence for junctions {first}/{second} "
                                 f"in {RETRIES} draws")
            out.append(dict(seq=seq, partner=reverse_complement(seq), sugarA=sugars, sugarB="R",
                            junction5=first, junction3=second))
    return out


def main():
    ensure_dirs()
    record = {}
    for modified, seed in SEEDS.items():
        rows = designs(modified, np.random.default_rng(seed))
        record[modified] = rows
        seen5 = Counter(r["junction5"] for r in rows)
        seen3 = Counter(r["junction3"] for r in rows)
        print(f"{modified}:R  {len(rows)} duplexes   "
              f"{modified}D|RR cells {len(seen5)}x{min(seen5.values())}-{max(seen5.values())}   "
              f"D{modified}|RR cells {len(seen3)}x{min(seen3.values())}-{max(seen3.values())}")
        print(f"        sugarA {rows[0]['sugarA']}   sugarB R")
        print(f"        example {rows[0]['seq']}  "
              f"junctions {rows[0]['junction5']} / {rows[0]['junction3']}")
    JUNCTION_RNA.write_text(json.dumps(record, indent=1))

    # append to the design the uniform arm started, so design.json describes the whole campaign
    rows = json.loads(DESIGN.read_text()) if DESIGN.exists() else []
    rows = [r for r in rows if r["family"] != "junction-rna"]
    for modified, entries in record.items():
        for i, r in enumerate(entries):
            rows.append(dict(name=f"JR{modified}_{i:03d}", family="junction-rna",
                             state=f"junction-{modified}:R", status="todo",
                             seq=r["seq"], partner=r["partner"],
                             sugarA=r["sugarA"], sugarB=r["sugarB"]))
    DESIGN.write_text(json.dumps(rows, indent=1))

    print()
    describe_design(rows)


if __name__ == "__main__":
    main()
