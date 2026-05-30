# AMBER MD pipeline for KLKB1 ASOs

This folder documents the AMBER molecular-dynamics pipeline currently used by
**Ariella Nouman** (Tamir lab) to simulate the ASO library targeting **KLKB1**.

The scripts themselves were authored by Ariella with the SLURM scaffolding
contributed by **Din**; they live on the TAU SLURM cluster at
`/tamir2/nouman/amber/`. They are mirrored here verbatim under
[`upstream_scripts/`](upstream_scripts/) so this repo can be the source of
truth for documentation and for proposed improvements. See
[`ATTRIBUTION.md`](ATTRIBUTION.md) for credit and licensing notes.

> **Status (2026-05-30):** first 5-PDB batch (`KLKB1_K1`) submitted from this
> repo on `power-general-shared-pool`. **3 of 5 jobs died at min1** with an
> OpenMPI core-binding error (`#processes: 32` vs SLURM cpuset on certain
> nodes). The 3 failures were retried as `KLKB1_K1_retry1` with `--bind-to
> none` added to every `mpirun` call in a v2 amber_pipe — see
> [`IMPROVEMENTS.md`](IMPROVEMENTS.md) §0 for the post-mortem.

---

## What the pipeline does

For each input PDB (a 20 nt ASO:RNA-or-DNA duplex, ~40 residues, pre-built
upstream), run the canonical AMBER setup → minimize → heat → equilibrate →
short production-MD chain:

1. **`tleap`** — load DNA.OL21 + RNA.OL3 + OPC water + the **modXNA** parameter
   sets for phosphorothioate (PS\*) and 2'-modified (FM\*, MC5, PC5, PM\*) chemistries;
   load the PDB; neutralize with Na⁺; solvate in an octahedral OPC box (+10 Å buffer);
   add ~150 mM NaCl; emit `${base}.prmtop` and `${base}.rst7`.
2. **Min1** (1000 cyc, 500 ncyc, restraints on residues 1-40 at 500 kcal/mol/Å²) —
   relax solvent + ions.
3. **Min2** (2500 cyc, 1000 ncyc, no restraints) — relax the full system.
4. **Heat** (50 000 steps × 2 fs = 100 ps; restraint 25 kcal/mol/Å² on ASO;
   linear ramp 100 → 310 K over the first 10 ps then hold).
5. **Equilibrate** (25 000 steps × 2 fs = 50 ps; NPT, 310 K, 1 atm;
   ASO restraint relaxed to 0.5 kcal/mol/Å²).
6. **Production MD** (1 000 000 steps × 2 fs = **2 ns**; NPT, 310 K, 1 atm,
   Langevin γ = 1.0 ps⁻¹, SHAKE on H bonds, frames every 1000 steps → 1000 frames).

Force-field / chemistry references (all loaded by `tleap`):
- `leaprc.DNA.OL21`, `leaprc.RNA.OL3`, `leaprc.water.opc` — standard Amber FFs.
- `/tamir2/nouman/amber/modXNA-main/dat/frcmod.modxna` and the per-residue
  `.lib` files (PSA/PSC/PSG/PST = phosphorothioate; FM\* = fluoro-modified;
  MC5, PC5, PM\* = 5'-cap / 2'-modified variants). These come from the
  **modXNA** project by **Erik Hartman et al.** (see
  [`ATTRIBUTION.md`](ATTRIBUTION.md)).

## Script layout (3 layers of bash)

```
jobs_req_new_Server_KLKB1[_P|_GPU|_GPU2].sh   # entry point — sbatch one job per PDB
        ↓ (each SLURM job runs)
amber_pipe_new_KLKB1.sh <pdb> <threads> <nodes>
        ↓ (delegates topology setup)
create_files_new_KLKB1.sh <pdb>               # tleap.in + 5 sander .in files
        ↓
tleap → sander.MPI ×5 (min1, min2, heat, equil, MD)
```

The four `jobs_req_...` variants only differ in their SLURM header:

| Script | Partition | QoS | Compute |
|---|---|---|---|
| `jobs_req_new_Server_KLKB1.sh` | (caller arg) | `owner` | CPU, `sander.MPI` |
| `jobs_req_new_Server_KLKB1_P.sh` | (caller arg) | `public` | CPU, `sander.MPI` |
| `jobs_req_new_Server_KLKB1_GPU.sh` | (caller arg) | `public` | 1× H100, `sander.MPI`* |
| `jobs_req_new_Server_KLKB1_GPU2.sh` | (caller arg) | `public` | 1× A100, `sander.MPI`* |

\* The GPU variants reserve a GPU but `amber_pipe_new_KLKB1.sh` defaults to
`USE_GPU=0` → it still launches `sander.MPI` on CPU. See
[`IMPROVEMENTS.md`](IMPROVEMENTS.md) §"GPU path is a no-op".

All four bake `name=$(whoami)` into the heredoc, so the SLURM job activates
`/tamir2/<runner>/miniconda3/envs/amber25` and runs from `/tamir2/<runner>/amber/`.
**This means a non-`nouman` runner needs a local copy of the 3 scripts** —
see [Running as a non-`nouman` user](#running-as-a-non-nouman-user).

## Running a batch

### As `nouman` (Ariella's pattern)

Inside `/tamir2/nouman/amber/`:

1. Create a small sibling folder, e.g. `KLKB1_A1/` (= 5 PDBs starting with `A`).
2. **Move** (not copy) PDBs from `KLKB1_run2/` into it. The move is the
   coordination primitive that prevents double-running when multiple people
   submit batches from the same source pool.
3. Submit:

   ```bash
   ./jobs_req_new_Server_KLKB1_P.sh KLKB1_A1 32 1 power-general-shared-pool
   #                                  ^batch  ^threads ^nodes  ^partition
   ```

   This sbatches one job per PDB. Each job requests `--mem=250G`, 32 cores ×
   1 node, `tamirtul-users_v2` account, up to 2 days.

### Running as a non-`nouman` user

If you can read but not write under `/tamir2/nouman/amber/`, copy the three
scripts into your own `/tamir2/$USER/amber/`, **copy** (not move) the PDBs into
your own batch folder, and patch the conda-env paths inside the `jobs_req_...`
heredoc to point at Ariella's env (which you can still read):

```bash
sed -i \
  -e 's|/tamir2/\${name}/miniconda3/condabin|/tamir2/nouman/miniconda3/condabin|g' \
  -e 's|/tamir2/\${name}/miniconda3/etc/profile.d/conda.sh|/tamir2/nouman/miniconda3/etc/profile.d/conda.sh|g' \
  -e 's|/tamir2/\${name}/miniconda3/envs/amber25|/tamir2/nouman/miniconda3/envs/amber25|g' \
  jobs_req_new_Server_KLKB1_P.sh
```

The `cd /tamir2/${name}/amber` line is left alone so amber_pipe + create_files
resolve to your local copies. The modXNA `.lib` and `frcmod` paths are
hardcoded to `/tamir2/nouman/amber/modXNA-main/` — these are world-readable so
they work as-is.

**When you copy (instead of move) from the shared pool, you MUST tell Ariella
which PDB names you took** so she doesn't queue the same names. The 2026-05-30
batch (`KLKB1_K1`) took:

- `ACGGTCTTCAAGCTGTTCTA.pdb`
- `ACTATAACAGTATCACTGTC.pdb`
- `ACTCAGGTTGTAAAAATTGC.pdb`
- `ACTGTCCTATATCACTCTAC.pdb`
- `ACTGTCCTATATCACTGTAC.pdb`

## Outputs (per ASO)

`create_files_new_KLKB1.sh` makes a subdirectory named after the PDB basename
and writes all generated files there. After the full pipeline finishes:

```
KLKB1_K1/ACGGTCTTCAAGCTGTTCTA/
├── ACGGTCTTCAAGCTGTTCTA_tleap.in
├── ACGGTCTTCAAGCTGTTCTA_min1.in / _min2.in / _heat.in / _equilibrate.in / _md.in
├── ACGGTCTTCAAGCTGTTCTA.prmtop / .rst7        # topology + initial coords
├── ACGGTCTTCAAGCTGTTCTA_min1.out / .ncrst     # minimization restart files
├── ACGGTCTTCAAGCTGTTCTA_min2.out / .ncrst
├── ACGGTCTTCAAGCTGTTCTA_heat.out / .ncrst / .nc
├── ACGGTCTTCAAGCTGTTCTA_equilibrate.out / .ncrst / .nc
└── ACGGTCTTCAAGCTGTTCTA_md.out / .rst7 / .nc  # production trajectory
```

The presence of `*_md.rst7` is currently the best (only) signal that the
pipeline completed end-to-end. There is no explicit `.complete` flag — see
[`IMPROVEMENTS.md`](IMPROVEMENTS.md) §"No success/failure marker".

## Monitoring & failure recovery (current state)

- `squeue -u $USER` shows queued/running jobs.
- Per-job stdout/stderr land in `slurm-<JOBID>.out` at the directory you
  invoked `jobs_req_...` from (NOT inside the batch folder).
- There is currently no automatic retry, no completion ledger, and no
  failed-job sweeper. Recovery is manual: read the `slurm-*.out`, identify the
  stage that died, and re-run.

See [`IMPROVEMENTS.md`](IMPROVEMENTS.md) for proposed automation around this.

## Coordination protocol (when running alongside Ariella)

From Ariella's own notes (paraphrased):

> Don't run aggressively over `KLKB1_run2/` as a whole. Carve off a sub-batch
> into a sibling folder (`KLKB1_A1`, `KLKB1_A2`, …), **MOVE** the PDBs out of
> `KLKB1_run2/`, and **don't leave duplicates** in the original pool. Anything
> still in `KLKB1_run2/` is "up for grabs". When a sub-batch falls over, move
> the failed PDBs into their own folder for retry. At the end, the results are
> aggregated.

If you cannot write to `KLKB1_run2/` (so cannot `mv`) you must **copy** out
and **manually tell Ariella the filenames you took** to preserve the
no-double-run invariant.
