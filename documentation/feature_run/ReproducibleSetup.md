# Reproducible setup

Path-relative setup for building the TAUSO environment reproducibly and running
the feature + model pipeline. **Everything lives under one base directory you
choose.** The commands assume nothing about how they're launched — run them
directly, or wrap each step with whatever job scheduler you use. No paths,
hostnames, or schedulers are hard-coded.

## Config

Pick one base directory; everything hangs off it:

```bash
export BASE=/path/to/base        # choose this once
```

Derived locations (used throughout):

| under `$BASE` | purpose |
|---------------|---------|
| `$BASE/TAUSO`        | the repo clone |
| `$BASE/bin`          | the micromamba binary |
| `$BASE/micromamba`   | env + package store (`MAMBA_ROOT_PREFIX`) |
| `$BASE/.tauso_data`  | downloaded data (`TAUSO_DATA_DIR`) |

Reproducibility anchors: the **explicit lock** (`conda-linux-64-dev.lock`) + a
**pinned git commit**. Record both with any run.

---

## Phase 0 — micromamba (isolated)

```bash
cd "$BASE"
curl -L -o micromamba \
  https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
chmod +x micromamba
mkdir -p "$BASE/bin"
mv micromamba "$BASE/bin/micromamba"

export PATH="$BASE/bin:$PATH"
export MAMBA_ROOT_PREFIX="$BASE/micromamba"

micromamba --version
micromamba info        # confirm root prefix points at $BASE/micromamba
```

For convenience you can add these two exports to your shell rc:

```bash
export PATH="$BASE/bin:$PATH"
export MAMBA_ROOT_PREFIX="$BASE/micromamba"
```

Do **not** run `micromamba shell init` — it injects an auto-activate hook.
`MAMBA_ROOT_PREFIX` only affects micromamba, so any other conda/envs you have
stay untouched.

---

## Phase 1 — clone, pinned, with submodules

```bash
cd "$BASE"
git clone --recurse-submodules https://github.com/RedPenguin100/TAUSO.git TAUSO
cd "$BASE/TAUSO"
git checkout <COMMIT>            # pin to a known commit and record it
git submodule update --init --recursive

git rev-parse HEAD
git submodule status             # both submodules present, no leading -/+
```

Submodules (public, HTTPS): `external/risearch`,
`notebooks/competitors/OligoAI/OligoAI`. No git-LFS.

---

## Phase 2 — env from explicit lock + editable install

Long commands are kept in scripts so they're trivial to launch any way you like.
[`build_env.sh`](./build_env.sh) **self-locates** — it derives the base as the
directory above the clone — and honours `MAMBA_ROOT_PREFIX` / `ENV` / `MM`
overrides:

```bash
bash "$BASE/TAUSO/documentation/feature_run/build_env.sh"
```

It creates the env `tauso_repro` from `conda-linux-64-dev.lock` (only if
missing), then editable-installs `tauso` (`uv` if present, else `pip`), then
verifies.

- Uses the `-dev` lock (includes `uv` + test/dev tooling — same env CI builds).
- Needs outbound internet (downloads packages).

---

## Phase 3 — data setup

```bash
export TAUSO_DATA_DIR="$BASE/.tauso_data"
mkdir -p "$TAUSO_DATA_DIR"
cd "$BASE/TAUSO"

# genome + bowtie + omics + raccess  (bowtie is the slow step — tune threads/mem)
micromamba run -n tauso_repro tauso setup-all -t 32 --mem-per-thread 2048

# cell context (cohort + expression + CAI weights + tGCN)
micromamba run -n tauso_repro tauso build-cell-context

# process + index the oligo data
micromamba run -n tauso_repro python -u -m notebooks.data.OligoAI.assign_canonical_gene
micromamba run -n tauso_repro python -u -m notebooks.utils.data
```

- `setup-all` is idempotent (hash-checks, skips existing; `--force` to rebuild).
- To isolate just the long pole:
  `micromamba run -n tauso_repro tauso setup-bowtie -t 32 --mem-per-thread 2048`
- Downloads genome / DepMap / bowtie sources — needs internet.

---

## Phase 4 — calculate all features

```bash
cd "$BASE/TAUSO"
micromamba run -n tauso_repro \
  python -u -m notebooks.features.calculate_features --dataset oligo --cpus 64 --overwrite
```

- `--overwrite` forces a fresh recompute; drop it to only fill gaps / resume.
- Heavy steps (`off_target_*`, `on_target_hybridization`, `structure`, `mfe`)
  can be split with `--step <name>` or `--partition K N` + `--merge N`.

---

## Phase 5 — models

```bash
cd "$BASE/TAUSO"
# example: gain + rank:ndcg
micromamba run -n tauso_repro \
  python -u -m notebooks.models.train_model --loss ndcg --split cohort --importance gain --cpus 64
```

---

Every command above is launch-agnostic and expressed relative to `$BASE`. Wrap
the long-running steps with your scheduler as you see fit; nothing here assumes
how the command is dispatched.
