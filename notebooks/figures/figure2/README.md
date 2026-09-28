# Figure 2 — TAUSO against the published tools

One script per panel. Each reads only what is in `data/`, computes its own numbers and writes its
own PNG, so a panel can be read, rerun and reviewed on its own.

```
build_data.py                    assembles data/ from the deployed model and the feature store
data/scores.parquet              the comparison subset: labels + every method's score
data/screening_effort.csv        per screen, the random-screening fold each method saves
data/parsimony.csv               the top-K feature sweep
panel_a_within_experiment.py     ranking inside one screen
panel_b_within_cohort.py         ranking across a gene x cell-line cohort
panel_c_shortlist.py             what each method's top 5% actually delivers
panel_d_by_chemistry.py          the same ranking, split by gapmer architecture
panel_e_parsimony.py             how much of the ranking survives dropping features
panel_f_screening_effort.py      wet-lab work saved against a random screen
out/                             rendered panels (not tracked)
```

## Redrawing

```bash
tauso setup-features --include-competition     # once: the competitor scores
python build_data.py                           # refresh data/ from the current model
python panel_a_within_experiment.py            # any panel, in any order
```

`build_data.py` is the only step that touches the model, the feature store or the training data.
Everything the panels draw is in `data/`, which is tracked: redrawing from a different model shows
up as a change to those files rather than silently under the plotting code.

## The comparison subset

Held-out test split, restricted to the two canonical gapmer architectures (5-10-5 2'-MOE and
3-10-3 cEt) and to rows every tool scored, so no method is credited or penalised for the rows it
happens to cover: 21,041 ASOs over 274 screens and 30 gene x cell-line cohorts.

A competitor score stored the other way round — larger meaning less knockdown — is flipped by
`build_data.py`, and only when it is anti-correlated in both groupings, which separates a reversed
convention from a tool that simply ranks badly in one of them. `OligoWalk·Tm` and `miRanda` are
flipped on the current data.

## `data/parsimony.csv`

Produced by a separate sweep that retrains on the top K gain-ranked features, which takes hours on
a GPU and is not run from this folder. The file shipped here is from the 485-feature model and
therefore does not match the other panels, which are drawn from the 679-feature deployed model;
panel e is stale until the sweep is rerun.
