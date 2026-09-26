# Open items

## Settle the equilibration cut-off across the whole campaign

`11_extract.py` discards the first 500 frames of each production run (`EQUILIBRATION_FRAMES`),
which is 1 ns while `prod.in` saves a frame every 2 ps. `check_frame_interval()` now enforces
that relationship, so the constant cannot quietly come to mean a different length of time.

Whether 1 ns is the right amount to discard has been checked twice, both times on a sample
too small to close the question.

What the checks showed, splitting production into half-nanosecond windows and comparing window
means against windows 5-8, on rep1 of RD_000, RM_000, RE_000, RR_000 and JRE_000:

- RR_000 and RM_000 are settled before the cut-off on every observable.
- RD_000 settles by 1.5 ns: its first kept window sits 3 sd out on hX and hRise, the next 0.1.
- RE_000 (hRise, hIncl) and JRE_000 (hX) are still moving into the second nanosecond.

The sizes are small. JRE_000 hX moves 0.5 A between the first kept window and the rest against
a dinucleotide spread of 4.8 A; RE_000 hRise moves 0.10 A against a spread of 0.6 A. So the
cut-off is not distorting the tables, but the tail of the discard window is not clean either,
and the drift has a consistent direction rather than averaging out.

To settle it:

- sample across all chemistries and both arms, and across replicates, not one run each
- as a script in `src/`, so the number is reproducible
- report it as "comparing observable means over successive windows", NOT "block analysis" --
  block averaging is the Flyvbjerg-Petersen error estimator and means something else

If the cut-off moves, `EQUILIBRATION_FRAMES` changes, all 1,920 runs need re-extracting, the
tables regenerate, and the frame count in the supplementary text (currently 1,500) changes.
