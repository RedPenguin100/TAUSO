# Open items

## Redo the equilibration cut-off properly

The 1 ns discard in `11_extract.py` (`EQUILIBRATION_FRAMES = 500`) rests on an ad-hoc test run
while the campaign was still going: 13 trajectories from 3 duplexes, all cEt junction systems
(`JRE_000/001/002`), by a throwaway script in a scratch directory.

What it showed: splitting each trajectory into four 1 ns windows, window 1 sat 2-3x the
between-window scatter away from windows 2-4, in the direction of the A-form starting geometry,
with no residual drift after it.

What to do instead, now that all 1,920 runs are extracted:

- sample across all chemistries and both arms, not one junction state
- enough runs for the estimate to be stable
- as a script in `src/`, not scratch, so the number in `11_extract.py` is reproducible
- report it as "comparing observable means over successive 1 ns windows", NOT "block analysis" --
  block averaging is the Flyvbjerg-Petersen error estimator and means something else

Consequence if the cut-off changes: `EQUILIBRATION_FRAMES` in `11_extract.py`, the frame count in
the supplementary text (currently 1,500), and every table in `~/md100` would need regenerating.
