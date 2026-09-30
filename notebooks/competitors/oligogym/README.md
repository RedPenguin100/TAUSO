# OligoGym baselines on the TAUSO split

[OligoGym](https://github.com/Roche/OligoGym) (Roche, Apache-2.0) is a benchmark suite for
oligonucleotide property prediction: curated datasets, featurizers that read HELM notation, and
classical + deep model wrappers. Its featurizers see only the oligo itself — sequence and
chemistry — with no genome, no accessibility and no off-target search.

That makes it the bar TAUSO's genome-derived features have to clear. This folder fits OligoGym's
classical models on the **train+val rows of TAUSO's frozen OligoAI split** and scores them on the
**same held-out test rows**, with the **same metric suite** (`notebooks/models/evaluate.py`), so
the numbers sit in one table next to the shipped `tauso_score_v1`.

## Files

| file | what it does |
|---|---|
| `helm_bridge.py` | TAUSO's `aso_sequence` / `chemical_pattern` / `ps_pattern` → a HELM string |
| `run_baselines.py` | fits each OligoGym featurizer × model on train+val, scores the test rows |
| `score_tauso.py` | scores the shipped TAUSO booster on those same test rows (trains nothing) |
| `score_competitors.py` | fetches the competition shards and scores all of them on those same rows |
| `compare.py` | merges both result files into one ranked table |
| `results/` | the generated `.txt` tables and `.json` records |

## Running it

Two environments, because OligoGym pins `scikit-learn==1.4.0` and pulls torch, lightning,
catboost, torch_geometric and tabpfn (`oligogym/models.py` imports all of them at module top,
even for the classical models):

```bash
# OligoGym side -- needs oligogym, not tauso
conda run -n oligogym_env python notebooks/competitors/oligogym/run_baselines.py

# reference rows -- need TAUSO_DATA_DIR (score_tauso.py also needs the tauso package)
conda run -n tauso_claude_ws python notebooks/competitors/oligogym/score_tauso.py
conda run -n tauso_claude_ws python notebooks/competitors/oligogym/score_competitors.py

conda run -n oligogym_env python notebooks/competitors/oligogym/compare.py
```

Each script puts the repo root and `src/` on `sys.path`, so the data path, the cohort definition,
the metric code and the column names all come from the repo rather than being restated here:

| taken from | what for |
|---|---|
| `notebooks.consts.OLIGO_CSV_PROCESSED_AVERAGED` | the dataset path |
| `notebooks.preprocessing.assign_cohort` | the `cohort_id` definition |
| `notebooks.models.evaluate.evaluate` | the metric suite |
| `tauso.data.consts` | every column name |

`notebooks/models/common.py` is the one that cannot be imported here — it reaches
`load_and_validate_final_data`, which pulls the whole feature pipeline (pyarrow, gffutils,
pyranges). So `METRICS`, the frozen `split` and the `reg`/`clean_exp` targets are still restated
in `run_baselines.py`, marked where they are. Installing the feature-pipeline dependencies into
`oligogym_env` would close that last gap.

`run_baselines.py` needs no feature store and no booster: everything it reads is in the processed
CSV. `oligogym_env` additionally needs `numba` and `pyarrow` for the two repo imports above.

## The HELM bridge

TAUSO stores chemistry as three parallel strings; HELM writes it as a dotted list of
`sugar(base)phosphate` monomers with no phosphate on the last residue:

```
RNA1{[cEt](G)[sp].[cEt](C)[sp].d(T)[sp].[cEt](A)}$$$$V2.0
```

The monomer vocabulary in `helm_bridge.py` was read off the HELM strings in OligoGym's own
datasets (`oligogym/resources/pkg_dataset/`). TAUSO's data uses only `M`/`C`/`d` sugars and
`*`/`d` linkages, which map to `[moe]`/`[cEt]`/`d` and `[sp]`/`p`.

Cross-check: OligoGym ships `Hwang_2024_1`, its own curation of the ASOptimizer collection, which
overlaps this dataset. On the 11,577 sequences present in both, **93.9% of the HELM strings this
bridge builds are byte-identical to OligoGym's own**. The remainder are gap-pattern disagreements
between the two curations of the same source ASOs, not encoding errors — TAUSO's chemistry
annotation is the reference here, since it is TAUSO's split being benchmarked.

## Result

Full table in `results/comparison.txt`. Train+val 121,571 rows, test 21,682. Top of the ranking:

| model | sees | exp_med | gxc_med | top5 | globP |
|---|---|---|---|---|---|
| TAUSO `tauso_score_v1` | 679 genome-derived features | **0.611** | **0.548** | **76.0** | **0.595** |
| OligoAI | sequence + transcript context | 0.435 | 0.320 | 67.5 | 0.461 |
| `kmers/xgboost/clean_exp` | sequence + chemistry | 0.345 | 0.330 | 62.5 | 0.346 |
| `OW_Duplex` (best non-OligoAI competitor) | OligoWalk duplex energy | 0.198 | 0.009 | 45.3 | -0.079 |
| `sfold_accessibility` | Sfold accessibility | 0.112 | 0.098 | 50.1 | 0.105 |
| `PFRED_PLS` | PFRED | 0.016 | 0.043 | 44.3 | -0.038 |
| `onehot/xgboost/clean_exp` | 0.329 | 0.287 | 61.5 | 0.318 |
| `thermo/xgboost/clean_exp` | 0.329 | 0.322 | 63.5 | 0.349 |
| `kmers/knn/reg` (worst-but-one) | 0.194 | 0.135 | 54.0 | 0.213 |

The best sequence-only baseline reaches 0.345 within-experiment Spearman against TAUSO's 0.611, a
gap of 0.266. OligoAI sits between them at 0.435, which is the point of including it: 0.345 is
where a bag-of-k-mers baseline belongs, not a broken harness. Roughly, 0.345 -> 0.435 is what
transcript context buys, and 0.435 -> 0.611 is what the accessibility, structure and off-target
features buy.

The competitor rows change how the baselines read: **every OligoGym configuration except the two
weakest beats every published tool other than OligoAI.** PFRED, OligoWalk, Sfold and miRanda all
sit at or below 0.20 exp_med, so a bag of k-mer counts is not a weak bar for this dataset -- it is
a stronger ranker than the established design tools. Only OligoAI, which is trained on this data's
own distribution, clears it.

Competitor scores are evaluated as they come, with no sign flipping: each tool's intended
direction is its own convention and guessing it would invent a result. A negative `exp_med` means
the column ranks opposite to inhibition; the magnitude is the signal it carries. `OW_Break_Target`
is constant over the test rows and is skipped, the same degeneracy check
`regen_common.validate_features` makes. Across the grid the model matters more than the featurizer: XGBoost and random
forest cluster at 0.31-0.35 whichever featurizer feeds them, while ridge lands at 0.15-0.27 and
k-nearest-neighbours at 0.14-0.20. Fitting the `clean_exp` target beats raw inhibition almost
everywhere, by roughly 0.01-0.04.

## What the comparison is and is not

- Both sides are scored on identical test rows by identical code, so the ranking is fair.
- The baselines are **sequence + chemistry only**, while TAUSO also has dose, cell line and
  transfection method. That sounds like an unfair edge but is not one for the headline metric:
  `volume_nm`, `cell_line` and `canonical_gene_name` are constant within every one of the 1,795
  experiments, and `exp_med` ranks within an experiment, so a feature that does not vary inside
  the group cannot change the within-group ranking. Those covariates move `globP`, and `gxc_*`
  a little, and `exp_med` not at all.
- Two targets are fitted, mirroring `notebooks/models/common.VARIANTS`: `reg` (raw
  `inhibition_percent`) and `clean_exp` (deviation from each experiment's mean, which is what
  TAUSO ships).
