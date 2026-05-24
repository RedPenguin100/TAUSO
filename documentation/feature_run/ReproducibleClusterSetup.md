# Reproducible cluster setup (powerslurm)

End-to-end, "nothing-polluted" setup for running the TAUSO feature + model
pipeline reproducibly on the SLURM cluster. Builds an isolated micromamba env
from the **explicit lock** (exact pinned packages, no solve) at a pinned git
commit, then runs every step via an explicit `micromamba run -n <env>` so it
never depends on the login PATH or `tauso_run`'s internal activation.

## Config (paths used below)

| var | value |
|-----|-------|
| work area      | `/tamir2/kovaliov` |
| micromamba bin | `/tamir2/kovaliov/bin/micromamba` |
| root prefix    | `/tamir2/kovaliov/micromamba`  (`MAMBA_ROOT_PREFIX`) |
| env name       | `tauso_repro` |
| repo checkout  | `/tamir2/kovaliov/TAUSO_fresh` |
| pinned commit  | `081be2dc7e4f329410080752c6b462a9605e3ea7` |

Reproducibility anchors: **the explicit lock** (`conda-linux-64-dev.lock`) +
**the pinned commit** above. Record both with any run.

---

## Phase 0 — install micromamba (isolated)  ✅

The `micro.mamba.run` tarball endpoint returned a corrupt stream on the login
node; the static binary from GitHub releases works:

```bash
cd /tamir2/kovaliov
curl -L -o micromamba \
  https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
chmod +x micromamba

mkdir -p /tamir2/kovaliov/bin
mv micromamba /tamir2/kovaliov/bin/micromamba

# verify
/tamir2/kovaliov/bin/micromamba --version
/tamir2/kovaliov/bin/micromamba -r /tamir2/kovaliov/micromamba info
```

### `.bashrc` convenience (safe — no shell init / no auto-activate)

```bash
# --- micromamba (isolated, for reproducible tauso env) ---
export PATH="/tamir2/kovaliov/bin:$PATH"
export MAMBA_ROOT_PREFIX="/tamir2/kovaliov/micromamba"
# optional shorthand:
# alias mmrun='micromamba run -n tauso_repro'
```

Do **not** run `micromamba shell init` — it injects an activation hook that can
auto-activate `base`. `MAMBA_ROOT_PREFIX` only affects micromamba, so your
`miniconda3` / `tauso_claude` are untouched. SLURM `srun`/`sbatch` inherit the
login shell env, so jobs see `MAMBA_ROOT_PREFIX` too.

---

## Phase 1 — fresh clone, pinned, with submodules  ✅

```bash
cd /tamir2/kovaliov
git clone --recurse-submodules https://github.com/RedPenguin100/TAUSO.git TAUSO_fresh
cd /tamir2/kovaliov/TAUSO_fresh
git checkout 081be2dc7e4f329410080752c6b462a9605e3ea7   # or stay on main tip
git submodule update --init --recursive

# verify
git rev-parse HEAD       # 081be2dc7e4f329410080752c6b462a9605e3ea7
git submodule status     # external/risearch + notebooks/competitors/OligoAI/OligoAI, no -/+ flags
```

Submodules (both public, HTTPS): `external/risearch`,
`notebooks/competitors/OligoAI/OligoAI`. No git-LFS.

---

## Phase 2 — env from explicit lock + editable install  ⏳ (run on a compute node)

Login node is slow for package extraction/linking, so build inside a job.

**Important `tauso_run` gotcha:** it does NOT pass complex commands (`&&` chains,
`|` pipes, `-n <env>` flags, quoted strings) through cleanly — it mangled an
inline one-liner, splitting `micromamba run -n` from its argument
(`tauso_repro: command not found`). **Fix: keep the logic in a script and run
`tauso_run bash <script>`** — two clean tokens, nothing to mangle. Use this
pattern for every later phase too.

The build script is committed at [`build_env.sh`](./build_env.sh) (so it's in
your checkout after Phase 1). Run it on a compute node:

```bash
tauso_run --cpu=8 --mem=32G \
  bash /tamir2/kovaliov/TAUSO_fresh/documentation/feature_run/build_env.sh
```

The script derives the repo root from its own location and honours `MM`,
`MAMBA_ROOT_PREFIX`, `ENV` overrides (defaults: the kovaliov paths + `tauso_repro`).
It creates the env from `conda-linux-64-dev.lock` only if missing, then editable-
installs `tauso` (`uv` if present, else `pip`), then verifies.

- Uses `conda-linux-64-dev.lock` (the `-dev` one — includes `uv` + test/dev tooling, same env CI builds).
- Output → the job's `slurm-<jobid>.out` (no `tee` pipe in the `tauso_run` line — same mangling reason).
- ⚠️ Needs outbound internet on the compute node (downloads packages). If the
  node is network-isolated, run the `create` on the login node instead, or use
  an internal mirror.
- **Status:** env `tauso_repro` created successfully (`micromamba create` finished);
  editable install via the script.

---

## Phase 3 — data setup (genome + bowtie + omics + oligo)  ⬜ pending

```bash
# genome + bowtie + omics + raccess (bowtie = slow step, custom thread params)
tauso_run --cpu=32 --mem=64G 'micromamba run -n tauso_repro tauso setup-all -t 32 --mem-per-thread 2048'

# cell context (cohort + expression + CAI weights + tGCN)
tauso_run --cpu=8 --mem=32G  'micromamba run -n tauso_repro tauso build-cell-context'

# process + index the oligo data
tauso_run --cpu=4 --mem=16G  'micromamba run -n tauso_repro python -u -m notebooks.data.OligoAI.assign_canonical_gene'
tauso_run --cpu=4 --mem=16G  'micromamba run -n tauso_repro python -u -m notebooks.utils.data'
```

- Bowtie alone, if you want to isolate the long pole:
  `tauso_run --cpu=32 --mem=64G 'micromamba run -n tauso_repro tauso setup-bowtie -t 32 --mem-per-thread 2048'`
- `setup-all` is idempotent (hash-checks, skips existing); add `--force` for a true rebuild.
- These download large files (genome, DepMap, bowtie sources) → same internet caveat.
- Ensure `TAUSO_DATA_DIR` is exported (target dir for all `setup-*` downloads), as CI does.

---

## Phase 4 — calculate all features  ⬜ pending

```bash
tauso_run --cpu=64 --mem=64G 'micromamba run -n tauso_repro python -u -m notebooks.features.calculate_features --dataset oligo --cpus 64 --overwrite'
```

- `--overwrite` forces a fresh recompute; drop it to only fill gaps / resume.
- Heavy steps (`off_target_*`, `on_target_hybridization`, `structure`, `mfe`)
  can be split with `--step <name>` or `--partition K N` + `--merge N`.

---

## Phase 5 — models  ⬜ pending

```bash
# example: gain + rank:ndcg (best ranking-quality config in local sweeps)
tauso_run --cpu=64 --mem=64G 'micromamba run -n tauso_repro python -u -m notebooks.models.train_model --loss ndcg --split cohort --importance gain --cpus 64'
```

---

## Monitoring jobs

```bash
my_status                          # squeue alias
tail -f /tamir2/kovaliov/TAUSO_fresh/env_build.log
# or the slurm-<jobid>.out file the job writes
```
