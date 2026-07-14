# MALAT1 2'-MOE ASOs ranked by OligoAI

Design 2'-MOE gapmer ASOs against **MALAT1** and rank them with the **OligoAI** inhibition model.
OligoAI scores a given ASO, so "design" is enumerate-then-score: TAUSO tiles the transcript into one
candidate per position, OligoAI scores each for predicted percent knockdown.

## Design
- 20-mer **5-10-5 MOE gapmer** (5 nt 2'-MOE / 10 nt DNA / 5 nt 2'-MOE), **full phosphorothioate**
- One candidate per transcript position → 8,810 candidates
- Scored under **Lipofection @ 100 nM** (the in-domain mode of OligoAI's Lipofection data)

Candidate generation uses TAUSO (`tauso.aso_generation.get_initial_data` on the MALAT1 transcript from
the genome cache), so it needs a TAUSO environment with `TAUSO_DATA_DIR` set.

## Run
From the repo root, in a TAUSO environment:
```bash
python -m notebooks.prediction.malat1_oligoai.design_malat1_oligoai   # --delivery / --dose to change context
```
All paths default relative to this file (no absolute paths) and are overridable by flag or env var:
- **OligoAI repo** — `--oligoai-repo` / `OLIGOAI_REPO`; defaults to the sibling `OligoAI-fork` checkout (holds `run_inference.py`).
- **Checkpoint** — `--ckpt` / `OLIGOAI_CKPT`; defaults to `checkpoints/OligoAI_11_09_25.ckpt` (gitignored, ~2.9 GB — put the trained OligoAI checkpoint there).
- **Inference env** — `--conda-env` / `OLIGOAI_ENV`; defaults to `oligo_5090_hybrid` (torch cu128 + flash-attn 2.8.4 for the RTX 5090).

If the OligoAI repo or checkpoint is missing, the script writes the OligoAI-format input and stops. The
OligoAI input and raw predictions go to a temp dir; only the ranked CSV is kept.

Off-target counts need the Bowtie index (`tauso setup-bowtie --genome GRCh38`). If it is missing the
count columns are filled with `<NA>` and scoring/ranking still complete — the build is not triggered.

## Output
`MALAT1_2moe_oligoai_ranked.csv` — 8,810 candidates sorted by `predicted_inhibition_percent`
(descending):
- `target_start` — 0-based transcript position of the site; `aso_sequence_5_to_3` — the ASO
- `predicted_inhibition_percent` — OligoAI predicted knockdown
- `perfect_matches` / `off_targets_1mm` / `off_targets_2mm` — genome-wide Bowtie off-target match
  counts at 0 / 1 / 2 mismatches (both strands, one pass), **excluding any hit inside the MALAT1
  locus** (the intended site and intragenic matches are not counted). So `perfect_matches == 0` means
  the ASO has no perfect match anywhere outside MALAT1 (96.6% of candidates); a non-zero value flags a
  real perfect off-target. This is a coarse sequence-specificity screen — see `tauso run-off-target`
  for the strand-aware, gene-level off-target detail.

## Notes
- OligoAI conditions only on sequence, sugar/backbone mods, ±50-nt context, dose, and delivery —
  **cell line, cell density, and volume are not model inputs**, so the ranking is cell-line-agnostic.
- MALAT1 is **in-domain** for OligoAI (2,915 MALAT1 rows in training): good for picking, not a blind benchmark.
