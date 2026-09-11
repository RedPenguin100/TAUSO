# DIAPH3 ordered panel — what was predicted before the experiment

Twelve ASOs ordered against DIAPH3 in NIH OVCAR-8, recorded before any of them were assayed.
The point of this file is the timestamp: the ranks and scores here are the model's predictions,
not a description written afterwards.

## What produced the numbers

| | |
|---|---|
| model | `tauso_score_v1_clean_exp_med.json`, md5 `ea59020955191a062ca35d3605665ca2` |
| deposit | Zenodo record 22688592, fetched by `tauso setup-model` |
| features | 670, `src/tauso/inference/model/tauso_score_v1.features.txt` |
| candidates | DIAPH3 pre-mRNA tiled into 20-mers, introns included: 498,327 rows, 489,337 unique sequences |
| tiling | `notebooks/prediction/design_DIAPH3.py`, deposited at Zenodo record 22693356 |
| context | OVCAR-8 (ACH-000696), lipofection, 100 nM |
| target | `clean_exp` -- inhibition centred within its screen, so +27.7 reads as 27.7 points above that screen's mean |

Ranks are over the 489,337 unique sequences, rank 1 = best predicted. Both backbones are
scored for every ASO, so `rank_full_ps` and `rank_5po` are both recorded whichever one is
being ordered.

## The panel

Eight predicted active, four lower down the list. Three are exonic (#4, #8, #10); the rest sit
inside introns, which follows from tiling the pre-mRNA -- 99% of the candidates are intronic.

Off-target counts are genome-wide Bowtie matches outside the DIAPH3 locus at 0, 1 and 2
mismatches. They are not part of the model's score and were run separately; two candidates
elsewhere in the tiling scored well yet matched over twenty thousand genomic sites, which is
why they are recorded here.

## Two things to read carefully afterwards

`#11` is labelled low in the order but the model does not predict it fails: rank 31,039 of
489,337 puts it at the 93.7th percentile, score +4.89. A knockdown there is the model being
right. `#12` at the 76th percentile is mildly positive too. The genuine low-end anchors are
`#9` (23rd percentile) and `#10` (36.5th).

Neither OVCAR-8 nor DIAPH3 appears in the training corpus, and intronic pre-mRNA targeting is
rare in it, so all twelve predictions are extrapolation rather than interpolation.
