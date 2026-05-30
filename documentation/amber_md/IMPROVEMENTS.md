# Proposed improvements to the AMBER MD pipeline

Ranked by **(impact × ease)**, with a candid take. Anything tagged *⚠ verify
empirically* is my guess from reading the scripts, not measured on the
cluster. None of these are urgent enough to block running the current
pipeline — they're investments that pay back if we end up running ≫ 184 ASOs
or want to publish the protocol.

## -1. Maintain a runtime `--exclude` list of wedged shared-pool nodes 🔴 *(empirically required — 2026-05-30)*

**Problem.** Some nodes on `power-general-shared-pool` are wedged in a way
that kills any job sent to them within 2 seconds, exit `0:53`
(`RaisedSignal:53 (Real-time signal 19)`), before `bash` even reads the
batch script. No `slurm-*.out`, no `.failed` marker, no MaxRSS. The 2026-05-30
run hit it on **compute-0-34, compute-0-73, compute-0-75, compute-0-80** —
every retry on the same node failed identically while
`compute-0-278` / `-338` / `-381` ran cleanly.

This matches the cluster skill's standing note:

> task wedges on specific shared nodes (e.g. compute-0-353/354) still
> happen — likely memory/IO contention. Mitigation: `--exclude=`
> problematic nodes when resubmitting.

**Fix.** On any resubmit, add
`#SBATCH --exclude=<comma-separated bad nodes>` to the heredoc. Keep
a runtime allowlist (or denylist) in a file the jobs_req reads, not
hardcoded in committed scripts — the bad-node set shifts over days.

**Not in the committed `robust/` scripts** for that reason. Apply inline
when needed.

## 0. Add `--bind-to none` to every `mpirun` 🔴 *(empirically required — 2026-05-30)*

**Problem.** On the first batch (KLKB1_K1, 5 PDBs on `power-general-shared-pool`),
**3 out of 5** jobs died at the very first `mpirun` (min1) with:

```
A request was made to bind that would require binding processes to more
cpus than are available in your allocation:
   Application: sander.MPI    #processes: 32    Binding policy: CORE
```

The SLURM allocation requested `--ntasks-per-node=32` and `mpirun -np 32`
asks OpenMPI to bind 32 ranks to 32 cores — but on some nodes in the
`power-general-shared-pool` partition, the cpuset given to the job has
fewer bindable cores than 32 (HT siblings, partition policy, or a shared
node with another job already pinning cores). The 2 jobs that landed on
"clean" nodes (compute-0-338, compute-0-381) ran fine.

**Fix.** Add `--bind-to none` to every `mpirun` in `amber_pipe_new_KLKB1.sh`:

```bash
mpirun --bind-to none -np ${total_procs} "$amber_exec".MPI -i ...
```

OpenMPI then spawns the ranks without trying to pin them to specific cores;
the Linux scheduler handles placement. Throughput cost is negligible for a
sander.MPI run of this size (it's bottlenecked by MPI comm and per-step
work, not by L1/L2 locality). This is the cheapest, most robust fix.

**Even better:** also drop `-np` to 8 or 16 (see §7) so the binding can't
hit the limit even on a constrained node.

**Status.** Mirrored in this repo's working copy under
`/tamir2/kovaliov/amber/amber_pipe_new_KLKB1_v2.sh` and used by
`jobs_req_new_Server_KLKB1_P_v2.sh` (not in `upstream_scripts/` — those are
verbatim from Ariella). Worth proposing upstream.

## 1. Add a success/failure marker per ASO 🟢 *(small change, big payoff)*

**Problem.** Today, the only signal that a job finished is the presence of
`${base}_md.rst7`. Finding failures across a batch means reading every
`slurm-*.out`.

**Fix.** Append to `amber_pipe_new_KLKB1.sh`:

```bash
# After the last successful mdrun step
touch "$output_path/.done"
trap 'touch "$output_path/.failed.${SLURM_JOB_ID:-noslurm}"' ERR
```

(Or, more carefully, write the failing stage name into `.failed`.)

A 5-line bash sweeper can then list ASOs missing `.done` and re-submit them.

## 2. Resume-from-stage logic 🟢

**Problem.** A node failure mid-MD wastes the whole 2-day allocation.

**Fix.** Before each `mpirun` call in `amber_pipe_new_KLKB1.sh`, skip the
stage if its `.ncrst` (or `.rst7` for the final MD) already exists *and* is
non-zero. This makes a re-submit on the same PDB pick up at the first
unfinished stage. The minimization + heat + equil stages are O(minutes) so
this is mostly about not wasting the production-MD slot.

```bash
[ -s "${base_name}_min1.ncrst" ] || mpirun -np $total_procs ... -i ..._min1.in ...
```

## 3. Right-size `--mem=250G` 🟢 *(empirically settled — 2026-05-30)*

**Problem.** A 20-mer DNA duplex (~50 atoms × 2) in a ~30 Å OPC box with 150
mM NaCl is on the order of 30 000–60 000 atoms total. `sander.MPI` typically
uses **1–4 GB** for boxes that size. `--mem=250G` is ~50–250× over-requested
and will block the job from landing on smaller nodes for no real reason.

**Fix.** Use `--mem=64G`.

**Empirical history (this batch):**
- `--mem=16G`: 4 of 5 jobs killed at allocation time by SLURM cgroup with
  `RaisedSignal:53 (Real-time signal 19)` — the `power-general-shared-pool`
  policy enforces a per-cpu memory floor that 16G/32-core trips.
- `--mem=64G`: empirically lands on shared-pool nodes without the cgroup
  kill, while still leaving 32-core sander.MPI massively over-provisioned.
- `--mem=250G` (upstream): works, but blocks low-memory nodes for no
  measurable benefit.

If a future benchmark with `sstat --format=MaxRSS` confirms actual MaxRSS
< 30 GB across all ASOs, 64G is the right number; if shared-pool policy
changes, revisit.

## 4. Fix the GPU path so it actually uses the GPU 🟡

**Problem.** `jobs_req_new_Server_KLKB1_GPU.sh` reserves an H100 (`--gres=gpu:H100:1`)
but the inner `amber_pipe_new_KLKB1.sh` defaults `USE_GPU=0`, so it runs
`sander.MPI` on CPU and the GPU sits idle. Even if `USE_GPU=1` is set, the
script invokes `pmemd.cuda.MPI` — `pmemd.cuda` is normally single-GPU and
doesn't use MPI. The MPI variant exists (multi-GPU) but rarely works without
a custom build.

**Fix.** In the GPU jobs_req heredoc, export `USE_GPU=1` and change
`"$amber_exec".MPI` → `"$amber_exec"` when `USE_GPU=1`. Confirm
`pmemd.cuda` exists in `/tamir2/nouman/miniconda3/envs/amber25/bin/`. A
single H100 will typically run 2 ns DNA-in-water in **minutes**, vs hours on
32 CPU cores. Massive speedup if the env supports it.

## 5. Tleap preflight on the login node 🟡

**Problem.** A malformed PDB (broken chain, missing atoms, wrong residue
names for modXNA) makes `tleap` fail ~30 s into the job, but the 32-core /
250 GB / 2-day allocation is still consumed for those 30 s. With 174 PDBs that
adds up to wasted slot-time and slower queue throughput.

**Fix.** Run `tleap -f <pdb>_tleap.in` on the login node before sbatching.
If it fails (non-zero exit or no `.prmtop` emitted), log the PDB to a
`bad_pdbs.txt` and skip the sbatch.

## 6. Add `#SBATCH --requeue` + signal trapping 🟡

**Problem.** If the scheduler preempts or the node fails, the job dies with
no record beyond `slurm-<id>.out` showing `CANCELLED`.

**Fix.** Add `#SBATCH --requeue` so SLURM auto-resubmits on node failure;
combined with the resume-from-stage logic (§2) this is self-healing for most
transient failures. Optionally trap SIGTERM to flush a checkpoint before
SLURM kills the job at time-limit.

## 7. CPU sweet-spot is probably much lower than 32 🟡 *(⚠ verify empirically)*

**Problem.** `sander.MPI` scaling for a 30-60 k atom box typically plateaus
around **4–8 cores** because of communication cost. 32 cores is likely giving
~1.2× the speed of 8 cores but locking 4× the slots.

**Fix.** Benchmark `-c 4 / 8 / 16 / 32` on one PDB and pick the knee. (The
cluster skill's own note "128 cores is the sweet spot" applies to TAUSO's MFE
workload — that's a totally different workload.)

## 8. MD length: 2 ns is short for any conformational claim 🟠 *(scope question)*

**Problem.** 2 ns of production MD is just barely past "system has finished
equilibrating." It's enough for short-range relaxation but not for any
binding/conformational analysis. Most ASO/duplex MD papers run 50–500 ns.

**Caveat.** This may be deliberate — if these runs are a *screening* signal
for downstream analysis (e.g. extracting a representative frame for docking,
or computing initial-frame RMSDs), 2 ns is defensible. If they're meant to
support a structural argument, it's too short.

**Suggested action.** Confirm with Ariella what the MD output is *used* for
downstream. If screening, document that and move on. If structural, extend
`nstlim` to ~50 000 000 (100 ns) and request more wall time.

## 9. Topology-build redundancy 🟠 *(⚠ verify)*

**Problem.** `create_files_new_KLKB1.sh` is re-invoked on every restart and
will rebuild the topology each time. If we add resume-from-stage (§2), we
should also skip the tleap step if `.prmtop` already exists. Otherwise a
re-submit silently rebuilds and any random-ion placement changes between runs
break trajectory continuity.

**Fix.** Inside `amber_pipe_new_KLKB1.sh`:

```bash
if [ ! -s "${base_name}.prmtop" ]; then
    ./create_files_new_KLKB1.sh "$pdb_file"
fi
```

This also makes the resume idempotent.

## 10. Per-ASO logging / structured run manifest 🟠

**Problem.** Discovering which jobid corresponds to which ASO requires
correlating SLURM accounting (`sacct`) with the slurm-out filename. There's
no per-batch manifest.

**Fix.** In `jobs_req_new_Server_KLKB1_P.sh`, after each `sbatch`, capture the
returned jobid and append a row to `<batch>/manifest.tsv`:

```
pdb_basename<TAB>jobid<TAB>submit_iso_ts<TAB>partition<TAB>account
```

Then a tiny aggregator at the end of the batch can join with `sacct` to get
state, wall time, MaxRSS, and the failure reason.

## 11. Move from shell to a thin Python CLI 🔴 *(nice-to-have; big refactor)*

**Problem.** Three layers of bash that pass arguments through heredocs and
shell-quote escapes are fragile. Variable expansion is split across two
shells (login + SLURM), which has already caused bugs (the GPU `USE_GPU`
default-0 → no-op, the `$(whoami)` bake-in).

**Fix.** Wrap the whole flow in a `tauso.amber` Python module:

```bash
python -m tauso.amber.run_batch KLKB1_K1 \
    --partition power-general-shared-pool \
    --threads 8 --mem 8G \
    --conda-env /tamir2/nouman/miniconda3/envs/amber25 \
    --modxna /tamir2/nouman/amber/modXNA-main \
    --use-gpu auto
```

Internal generation of `.in` files, sbatch submission, resume detection,
manifest writing, and per-stage logging all live in Python with proper tests.
Cost: a couple of days. Payoff: a re-usable, version-controlled MD entry
point for the next target after KLKB1.

## 12. Capture/validate the `addions … 150 / 150` molarity 🔴 *(minor)*

**Problem.** `addionsrand dna1 Na+ 150 / Cl- 150` adds a fixed *count*, but
the actual molarity depends on the (variable) box volume per ASO. So the
nominal "150 mM NaCl" is system-dependent.

**Fix.** Compute the box volume after `solvateoct` and add ions to hit a
target *concentration* instead of a fixed count. (`packmol` or
`tleap`'s own ion-by-concentration trick.)

## 13. Documented force field versions 🔴 *(housekeeping)*

The pipeline uses `DNA.OL21` + `RNA.OL3` + OPC + modXNA. The exact versions
in Ariella's `amber25` env should be pinned in this README (or an
`environment.yaml`) so we can reproduce in five years.

## Summary table

| # | Improvement | Effort | Impact | Risk |
|---|---|---|---|---|
| -1 | `--exclude=<wedged nodes>` on shared pool | per-retry | **Required when nodes wedge** *(saw 4/5 jobs killed at 2 s)* | List rots — keep runtime-only |
| 0 | `--bind-to none` on mpirun | 2 min | **Required** *(60% job failure rate without it on the first batch)* | None |
| 1 | `.done`/`.failed` markers | 5 min | High | None |
| 2 | Resume-from-stage | 15 min | High | Low |
| 3 | Right-size memory | 10 min + 1 bench job | Medium-High | None |
| 4 | Fix GPU path | 30 min + GPU bench | High *(if it works)* | Medium |
| 5 | Tleap preflight | 15 min | Medium | None |
| 6 | `--requeue` + signal trap | 10 min | Medium | Low |
| 7 | Right-size CPU count | 30 min bench | Medium | None |
| 8 | Extend MD length | 0 code | Depends | Scope question |
| 9 | Idempotent topology build | 5 min | Medium *(only with §2)* | None |
| 10 | Per-batch manifest | 30 min | Medium | None |
| 11 | Python CLI rewrite | 2 days | High long-term | Migration cost |
| 12 | Ion-by-concentration | 1 hr | Low | None |
| 13 | Pin FF versions | 1 hr | Low (reproducibility) | None |

**Recommended near-term:** §1, §2, §6, §9 together (an afternoon) get us
self-healing batches. §3 + §7 together (one benchmark batch) cut the
per-ASO resource ask 4–8×, which means more parallel jobs land sooner.
§4 (GPU) is the single biggest speed win if `pmemd.cuda` is in the env.
