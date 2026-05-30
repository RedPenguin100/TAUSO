# robust/

Hardened drop-in replacements for the `upstream/` KLKB1 scripts.

- `amber_pipe_new_KLKB1_v3.sh` — replaces `amber_pipe_new_KLKB1.sh`.
- `jobs_req_new_Server_KLKB1_P_v3.sh` — replaces `jobs_req_new_Server_KLKB1_P.sh`.
- `create_files_new_KLKB1.sh` — verbatim copy of upstream, included so
  this folder is self-contained when deployed.

See `../IMPROVEMENTS.md` for which items are baked in (§0, §1, §2, §3, §5,
§6, §9).

## Deploy

    rsync -av notebooks/molecular_dynamics/robust/ tauso:/tamir2/$USER/amber/

The `_v3` suffix means these coexist safely with any upstream v1/v2 you
already have on the cluster.

## Run a batch

From `/tamir2/$USER/amber/`:

    ./jobs_req_new_Server_KLKB1_P_v3.sh KLKB1_K1 32 1 power-general-shared-pool

What changes vs. upstream:
1. Each PDB gets a tleap preflight on the login node. Bad PDBs land in
   `KLKB1_K1/bad_pdbs.txt` and never reach SLURM.
2. The SLURM job sees `.prmtop` already exists and skips topology rebuild.
3. Each stage is gated on its restart file (`*_min1.ncrst`,
   `*_min2.ncrst`, …). A re-submit picks up at the first unfinished stage.
4. Every `mpirun` has `--bind-to none`.
5. `--mem=64G`, `#SBATCH --requeue`.

## Find what failed / what's done

    find KLKB1_*/ -mindepth 2 -maxdepth 2 -name '.failed.*'
    find KLKB1_*/ -mindepth 2 -maxdepth 2 -name '.done'

## Retry

Re-run the same `jobs_req_v3` command on the batch directory. Done ASOs
are no-ops; partial ones resume at their first unfinished stage. No folder
gymnastics.

If wedged nodes hit (jobs exit `0:53` in 2 s with no `slurm-*.out`), add
`#SBATCH --exclude=<comma-separated node list>` to the heredoc in
`jobs_req_v3` for the resubmit. The bad-node list shifts day to day — see
`../IMPROVEMENTS.md` §-1.
