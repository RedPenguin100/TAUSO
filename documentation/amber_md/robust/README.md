# `robust/` — our hardened version of the pipeline

These are the scripts we actually run. They replace the upstream copies in
[`../upstream_scripts/`](../upstream_scripts/) and bake in the
robustness/efficiency improvements documented in
[`../IMPROVEMENTS.md`](../IMPROVEMENTS.md):

| Item | What it does | Where |
|---|---|---|
| §0 | `--bind-to none` on every `mpirun` | `amber_pipe_v3` |
| §1 | `.done` / `.failed.<jobid>.<stage>` markers | `amber_pipe_v3` |
| §2 | Resume-from-stage (skip a stage if its `.ncrst` exists) | `amber_pipe_v3` |
| §3 | `--mem=64G` (down from `250G`; 16G was empirically too low — shared-pool cgroup kill) | `jobs_req_v3` |
| §5 | Tleap preflight on login node, log bad PDBs, skip sbatch | `jobs_req_v3` |
| §6 | `#SBATCH --requeue` + SIGTERM/SIGINT trap | both |
| §9 | Skip `create_files` + `tleap` if `.prmtop` already exists | `amber_pipe_v3` |

Skipped (deliberately, for now): §4 (GPU path — needs separate verification
that `pmemd.cuda` is in the env) and §7 (cpu count benchmark — keeping the
old `-c 32` until we measure).

## Files

| File | Same as upstream? | Notes |
|---|---|---|
| `create_files_new_KLKB1.sh` | **Yes**, verbatim copy | Kept here so `robust/` is self-contained when deployed. |
| `amber_pipe_new_KLKB1_v3.sh` | No — see header comment | Drop-in replacement for `amber_pipe_new_KLKB1.sh`. |
| `jobs_req_new_Server_KLKB1_P_v3.sh` | No — see header comment | Drop-in replacement for `jobs_req_new_Server_KLKB1_P.sh`. |

## Deploy

Copy all three to your own `/tamir2/$USER/amber/` (alongside any in-flight
v1/v2 scripts; the `_v3` suffix means they coexist safely):

```bash
rsync -av documentation/amber_md/robust/ tauso:/tamir2/$USER/amber/
chmod +x /tamir2/$USER/amber/*.sh
```

## Run a batch

From `/tamir2/$USER/amber/`:

```bash
./jobs_req_new_Server_KLKB1_P_v3.sh KLKB1_K1 32 1 power-general-shared-pool
```

What changes vs. upstream when you run this:

1. **Login-node preflight** — for each PDB it generates the `.in` files and
   runs `tleap` on the login node. Takes ~30 s per PDB. If tleap fails the
   PDB is appended to `KLKB1_K1/bad_pdbs.txt` and no SLURM job is submitted.
2. **The SLURM job sees `.prmtop` already exists** (from preflight) and skips
   `create_files` + `tleap`, going straight to `min1`.
3. **Each stage is gated** — if a previous run wrote `min1.ncrst`, this run
   skips `min1` and starts at `min2`. Same for every stage through MD.
4. **On failure** — a `.failed.<jobid>.<stage>` marker lands in the ASO's
   subfolder. On success, all `.failed.*` are cleared and `.done` is written.
5. **On node failure or preemption** — SLURM auto-resubmits (`--requeue`).
   The resubmitted job sees the partial `.ncrst` files and resumes there.

## Find what failed

```bash
# All failed ASOs across all batches
find KLKB1_*/ -mindepth 2 -maxdepth 2 -name '.failed.*'

# All successful ASOs
find KLKB1_*/ -mindepth 2 -maxdepth 2 -name '.done'

# Failed at a specific stage
find KLKB1_*/ -mindepth 2 -maxdepth 2 -name '.failed.*.min1'
```

## Retry a failed batch

Just re-run `jobs_req_v3` on the same batch directory. The preflight will
skip already-built topologies, and `amber_pipe_v3` will skip already-completed
stages. ASOs whose `.done` already exists are no-ops; ASOs whose last `.ncrst`
is `heat.ncrst` will resume at `equilibrate`. No folder moves needed.

## What is **not** changed

- The force fields / `tleap` setup (in `create_files_new_KLKB1.sh`) — verbatim.
- The simulation parameters (`min1/min2/heat/equilibrate/md` `.in` files) —
  verbatim. Same 2 ns production MD.
- The `name=$(whoami)` heredoc convention — the SLURM job still `cd`s to
  `/tamir2/<runner>/amber/`. The conda env is hardcoded to Ariella's
  `/tamir2/nouman/miniconda3/envs/amber25` since that's where the install
  lives.
