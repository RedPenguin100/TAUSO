# Improvements

Status: ✅ done in `robust/` · ⏭ skipped intentionally · ◇ open / not yet done.

| # | Item | Status |
|---|------|--------|
| -1 | `#SBATCH --exclude=<wedged nodes>` on shared-pool retries | runtime only — list rots, not committed |
| 0  | `--bind-to none` on every `mpirun` | ✅ |
| 1  | `.done` / `.failed.<jobid>.<stage>` markers per ASO | ✅ |
| 2  | Resume-from-stage (skip a stage if its `.ncrst` exists) | ✅ |
| 3  | `--mem=64G` (was 250G) | ✅ |
| 4  | GPU path (`pmemd.cuda`) — upstream GPU variants reserve a GPU then run sander.MPI on CPU anyway | ⏭ — needs verifying `pmemd.cuda` is in env |
| 5  | Tleap preflight on the login node; bad PDBs → `bad_pdbs.txt`, skip sbatch | ✅ |
| 6  | `#SBATCH --requeue` + SIGTERM trap | ✅ |
| 7  | Right-size `-np` (probably 4–8 enough; `sander.MPI` plateaus early) | ⏭ — needs benchmark |
| 8  | 2 ns is short for conformational claims — extend to 50–500 ns | ◇ — scope question, ask Ariella what the trajectories feed |
| 9  | Skip `create_files` + `tleap` if `.prmtop` exists | ✅ |
| 10 | Per-batch manifest (pdb, jobid, partition, submit_ts) | ◇ |
| 11 | Rewrite the 3-layer bash as a thin `tauso.amber` Python CLI | ◇ |
| 12 | Ions by concentration, not fixed count | ◇ |
| 13 | Pin force-field versions in an `environment.yaml` | ◇ |

## Empirical learnings (2026-05-30 run)

**§-1 — Wedged shared-pool nodes.** Specific nodes on
`power-general-shared-pool` kill any job sent to them within 2 s with exit
`0:53` (`RaisedSignal:53 / Real-time signal 19`), before bash reads the
script. No `slurm-*.out`, no marker, no MaxRSS. Today: `compute-0-34, 73,
75, 80`. Mitigation: inline `#SBATCH --exclude=<list>` on retries. Don't
commit the node list — it shifts day to day.

**§0 — OpenMPI core-binding mismatch.** Without `--bind-to none`, 3 of 5
v1 jobs died at min1 with `Application: sander.MPI #processes: 32` because
SLURM's cpuset gave fewer bindable cores than 32 on some shared-pool
nodes. Adding `--bind-to none` removed the constraint and worked on every
node we landed on after.

**§3 — `--mem` floor on the shared pool.** `--mem=16G` triggered the same
cgroup-level `RaisedSignal:53` kill (separate cause from §-1). `--mem=64G`
worked. Still ~4× under upstream's 250G, still ~16× over sander.MPI's
actual MaxRSS for a 20-mer in OPC. If a real benchmark shows MaxRSS <
30 GB, 64G is the right number.
