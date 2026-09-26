# Open items

## Widen the sample behind the equilibration cut-off

`11_extract.py` discards the first 500 frames of each production run (`EQUILIBRATION_FRAMES`),
which is 1 ns while `prod.in` saves a frame every 2 ps. `check_frame_interval()` enforces that
relationship, so the constant cannot come to mean a different length of time.

The cut-off holds. Splitting production into half-nanosecond windows and comparing window means
against windows 5-8, on rep1 of RD_000, RM_000, RE_000, RR_000 and JRE_000:

- RR_000 and RM_000 are settled before the cut-off on every observable.
- RD_000 settles by 1.5 ns.
- RE_000 (hRise, hIncl) and JRE_000 (hX) take a little longer, by about a tenth of their
  dinucleotide spread: JRE_000 hX moves 0.5 A where the spread is 4.8 A, RE_000 hRise 0.10 A
  where the spread is 0.6 A.

That is small enough to leave the tables as they are. What is worth doing when convenient is
running the same comparison across replicates and more duplexes, as a script in `src/` rather
than by hand, so the number is reproducible from the repository.

Report it as "comparing observable means over successive windows", NOT "block analysis" --
block averaging is the Flyvbjerg-Petersen error estimator and means something else.

If the cut-off ever moves, `EQUILIBRATION_FRAMES` changes, all 1,920 runs need re-extracting,
the tables regenerate, and the frame count in the supplementary text (currently 1,500) changes.
