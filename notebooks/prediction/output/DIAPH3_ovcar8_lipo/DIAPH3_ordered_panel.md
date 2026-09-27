# DIAPH3 ordered panel — scored by the 679-feature model

The twelve ASOs ordered against DIAPH3 in NIH OVCAR-8, re-scored after the shape-junction and
md100 work. The ASOs are the ones physically ordered; the ranks and scores here are the current
model's, not the ones that selected them. The predictions that predate the assay are the
670-feature numbers, in this file's git history.

## What produced the numbers

| | |
|---|---|
| model | `tauso_score_v1_clean_exp_med.json`, md5 `0f0364451aebbaff98115a8367510694` |
| deposit | registered in `src/tauso/inference/model/models.json`; not yet on a public deposit |
| features | 679, `src/tauso/inference/model/tauso_score_v1.features.txt` |
| candidates | DIAPH3 pre-mRNA tiled into 20-mers, introns included: 498,327 rows, 489,337 unique sequences |
| tiling | `notebooks/prediction/design_DIAPH3.py`; the 670-feature tilings are at Zenodo record 22693356 |
| context | OVCAR-8 (ACH-000696), lipofection, 100 nM |
| target | `clean_exp` -- inhibition centred within its screen, so +25.9 reads as 25.9 points above that screen's mean |

Ranks are over the 489,337 unique sequences, rank 1 = best predicted, and each ASO carries the
rank and score for the backbone it was ordered with: `#1`-`#4` and `#9`-`#12` full PS,
`#5`-`#8` 5 PO. A backbone changes 9 of the 679 features, so a rank only means anything
within its own backbone.

## The panel

Eight predicted active, four lower down the list. Three are exonic (#4, #8, #10); the rest sit
inside introns, which follows from tiling the pre-mRNA -- 99% of the candidates are intronic.

Off-target counts are genome-wide Bowtie matches outside the DIAPH3 locus at 0, 1 and 2
mismatches. They are not part of the model's score, were run separately, and are carried over
unchanged: they depend on the sequence and the genome, not the model.

## Two things to read carefully afterwards

`#11` is labelled low in the order but the model does not predict it fails: rank 25,181 of
489,337 puts it at the 94.9th percentile, score +5.26. A knockdown there is the model being
right. `#12` at the 75.5th percentile is mildly negative. The genuine low-end anchors are
`#9` (24.3rd percentile) and `#10` (30.1st).

The strongest pick in each backbone is `#3`, first of 489,337 on full PS, and `#7`, third on
5 PO.

Neither OVCAR-8 nor DIAPH3 appears in the training corpus, and intronic pre-mRNA targeting is
rare in it, so all twelve predictions are extrapolation rather than interpolation.
