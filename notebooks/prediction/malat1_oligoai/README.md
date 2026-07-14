# MALAT1 2'-MOE ASOs ranked by OligoAI

Design 2'-MOE gapmer ASOs against **MALAT1** (nuclear lncRNA) and rank them with the
**OligoAI** inhibition model. OligoAI is a *scorer*, not a generator, so this enumerates one
candidate per transcript position and scores each for predicted percent knockdown.

## Design (fixed)
- **20-mer, 5-10-5 MOE gapmer**: 5 nt 2'-MOE wing / 10 nt DNA gap / 5 nt 2'-MOE wing (RNase-H mechanism)
- **Full phosphorothioate** backbone (19 PS linkages)
- One candidate per position → the MALAT1 sense transcript (8,829 nt) yields **8,810** candidates
- Scored under **Lipofection @ 100 nM** (both are OligoAI model inputs; 100 nM is the in-domain Lipofection mode)

## Inputs
- `MALAT1.fasta` — MALAT1 sense transcript. Reconstructed from the tiled design and validated
  base-for-base against OligoAI's own MALAT1 context windows (2,905/2,905 rows; 80/80 20-mers exact).
- OligoAI checkpoint (**not committed**, ~2.9 GB): `OligoAI_11_09_25.ckpt` (trained RiNALMo-giga).
- OligoAI repo (`run_inference.py`, model code) and its conda env.

## Run
```bash
python design_malat1_oligoai.py --ckpt /path/to/OligoAI_11_09_25.ckpt
# quick partial run: add --limit 50
# other delivery contexts: --delivery Gymnosis|Electroporation --dose <nM>
```
Inference runs inside the `oligo_5090_hybrid` conda env (torch cu128 + flash-attn 2.8.4, required
for the RTX 5090). The script tiles the FASTA, writes OligoAI-format input, calls `run_inference.py`,
and ranks by predicted knockdown.

## Output
- `MALAT1_2moe_oligoai_ranked.csv` — **the deliverable**: 8,810 candidates ranked by
  `predicted_inhibition_percent`, with `target_start` (0-based position on the transcript) and sequence.
- Intermediates (`*_input.csv`, `*.with_predictions.csv`) are regenerable and git-ignored.

## Caveats
- MALAT1 is **in-domain** for OligoAI (2,915 MALAT1 rows in its training data), so these scores are for
  *picking* candidates, not a blind-generalization benchmark.
- Delivery method and dose are model inputs — changing them shifts the ranking.
- OligoAI parameterizes MOE and DNA; this design uses only those, so it is squarely in its modeled chemistry.
