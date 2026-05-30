# molecular_dynamics

AMBER MD scripts for the KLKB1 ASO library on the TAU SLURM cluster.

- `upstream/` — Ariella's KLKB1-family scripts, verbatim from `/tamir2/nouman/amber/`.
- `robust/` — hardened drop-in replacements (resume-from-stage, success/failure
  markers, login-node tleap preflight, `--bind-to none`, right-sized memory,
  `--requeue`). This is what we run. See `robust/README.md`.
- `IMPROVEMENTS.md` — list of improvements with their status.

The flat files at this level (`amber_pipe.sh`, `create_files.sh`, `jobs_req.sh`,
`pdb_files/`) are an earlier sketch unrelated to the KLKB1 work.
