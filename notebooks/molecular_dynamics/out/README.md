# out/ -- MD output, not tracked

A local mirror of the finished runs the lookup tables were built from, so the tables can be
rebuilt without going back to the cluster.

    systems/<name>/                  384 designs: RD RM RE RR (uniform) + JRM JRE (junction)
    systems/<name>/rep{1..5}/        5 replicates each, 1,920 runs in all
    runs/rodu_rna/rep1/<name>/       the earlier 2 ns cross-chemistry set
    runs/junction_rna/rep{1,2,3}/<name>/

Each design directory holds `sys.prmtop` and `sys.rst7` from solvation; each replicate holds the
minimisation, heating and equilibration logs, `md_nowat.nc` (4 ns of production at a frame every
2 ps, water and ions stripped from the coordinates), and the `features.csv` that step 11 wrote
from it.

About 15 GB in all. Everything here is reproducible by re-running `src/` on the cluster, which is
why only the tables below are tracked:

    features/tables/{step,pair,pair_triplet,res}.csv   the per-dinucleotide lookup tables
    features/screen.csv                                which observables were kept, and why
    features/correlations.csv                          the pairwise correlations behind it
