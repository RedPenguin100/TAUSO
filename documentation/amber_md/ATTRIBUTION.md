# Attribution

## Scripts mirrored in `upstream_scripts/`

The six scripts under [`upstream_scripts/`](upstream_scripts/) were authored by
**Ariella Nouman** (Tamir Tuller lab, Tel Aviv University), with the SLURM
job-submission scaffolding contributed by **Din** (Tamir lab member, full name
TBD). They are mirrored here **verbatim** from
`/tamir2/nouman/amber/` on the TAU SLURM cluster as of 2026-05-30. No
modifications have been made to the upstream copies in this folder; any
patches needed to run them as a non-`nouman` user are documented in
[`README.md`](README.md) and not applied to the mirror.

If we ever modify these scripts in-place under `upstream_scripts/`, the
modification log belongs in a sibling `CHANGES.md` so the mirror stays
auditable against Ariella's working copy.

## Force-field libraries (modXNA)

`create_files_new_KLKB1.sh` loads parameter files from
`/tamir2/nouman/amber/modXNA-main/` (`frcmod.modxna` + per-residue `.lib`
files for PS\*, FM\*, MC5, PC5, PM\*). These are **not** in this repo — they
come from the **modXNA** project (Erik Hartman et al.). Upstream:
<https://github.com/ErikHartman/modXNA> (verify the exact URL — TODO). The
`README.md` and `LICENSE` files in `modXNA-main/` on the cluster carry the
upstream's own terms; we should not redistribute them here without checking
that license.

## Core software

- **AmberTools / Amber25** — used for `tleap`, `sander.MPI`, `pmemd.cuda`.
  Licensed via the Amber license (GPLv3 for AmberTools, separate for
  pmemd/sander). Installed under
  `/tamir2/nouman/miniconda3/envs/amber25/`.
- **Amber force fields** — `DNA.OL21` (Galindo-Murillo, Šponer, Otyepka,
  Jurečka, Cheatham), `RNA.OL3` (Zgarbová, Otyepka, Šponer, Mládek, Banáš,
  Cheatham, Jurečka), water `OPC` (Izadi, Anandakrishnan, Onufriev). These
  ship with AmberTools.

## Pending audit items (see project memory)

The [attribution-audit-todo](../../../../../home/michael/.claude/projects/-mnt-datadrive-workdir-aso-project-TAUSO-research/memory/attribution-audit-todo.md)
memory tracks all open attribution questions for the TAUSO repo. When the
modXNA upstream URL/license is confirmed, update that memory and this file.
