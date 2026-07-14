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
From the repo root (the script imports `notebooks.data.OligoAI`):
```bash
export OLIGOAI_REPO=/path/to/OligoAI      # holds run_inference.py + the rinalmo package
export OLIGOAI_CKPT=/path/to/OligoAI.ckpt # trained checkpoint (not committed, ~2.9 GB)
python -m notebooks.prediction.malat1_oligoai.design_malat1_oligoai  # --delivery / --dose to change context
```
Without `OLIGOAI_REPO`/`OLIGOAI_CKPT` the script writes the OligoAI-format input and stops. Inference
runs in the conda env named by `OLIGOAI_ENV` (default `oligo_5090_hybrid`). The OligoAI input and raw
predictions are written to a temp dir; only the ranked CSV is kept.

## Output
`MALAT1_2moe_oligoai_ranked.csv` — 8,810 candidates ranked by `predicted_inhibition_percent`, with
`target_start` (0-based transcript position) and the ASO sequence.

## Notes
- OligoAI conditions only on sequence, sugar/backbone mods, ±50-nt context, dose, and delivery —
  **cell line, cell density, and volume are not model inputs**, so the ranking is cell-line-agnostic.
- MALAT1 is **in-domain** for OligoAI (2,915 MALAT1 rows in training): good for picking, not a blind benchmark.
