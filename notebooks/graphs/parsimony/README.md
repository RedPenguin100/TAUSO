# Feature-count descent (679 to 10 features)

Starting from the 679-feature matrix, the lowest-gain features were dropped step by step:
679 down to 100 in steps of 5, then 100 down to 10 one feature at a time (207 feature counts in all).

- **CV** at each feature count: 15 fits (3 fold-columns x 5 gene-grouped folds); the reported number is the mean over the 15 fits, `_sd` columns are the sd.
- **TEST** at each feature count: trained on all train+val rows, scored on the frozen 21,682-row TEST split with `evaluate`, mean over 3 seeds.
- `colsample_bytree = max(0.17189541, min(1.0, 40/n))`; the other parameters are the MED `clean_exp` deploy config.
- OligoAI is scored on the same TEST rows with the same `evaluate`: exp_med 0.4349, gxc_med 0.3205.

## Files

| file | what it is |
|---|---|
| `results_by_n.csv` | one row per feature count: CV and TEST rmse, mae, exp_med, exp_mean, gxc_med, gxc_mean, and TEST top5, p5, p10, globP |
| `curve.csv`, `test_curve.csv` | the same numbers with their sds, per feature count, for CV and TEST |
| `drop_order.csv` | every feature, the last n it was in the model and the first n without it; features dropped in the same step share one row of n (`dropped_together_with` says how many others left in that step); 10 features are never dropped |
| `features_by_n.json` | the exact feature list at each feature count, as `{"679": [...], ..., "10": [...]}` |
| `other_models_test.csv` | TEST exp_med for the LOW/MED `clean_exp` and `regression` deploy configs on the same matrix and split, 3 seeds each |
| `oligoai_test_panel.json` | OligoAI's TEST scores |
| `make_tables.py` | rebuilds `results_by_n.csv` and `drop_order.csv` |
| `plot_vs_oligoai.py` | `tauso_vs_oligoai_parsimony.png`: TEST exp_med against feature count, with OligoAI and the MED configs as lines |
| `plot_curve.py` | `parsimony_curve.png`: CV (top) and TEST (bottom) for rmse, exp_med, exp_mean, gxc_med, with standard-error bars |

## Results

| n | CV rmse | CV exp_med | TEST rmse | TEST exp_med | TEST gxc_med | TEST top5 |
|---|---|---|---|---|---|---|
| 679 | 19.0921 | 0.5839 | 19.6803 | 0.6126 | 0.5524 | 75.75 |
| 644 | 19.0886 | 0.5849 | 19.6676 | 0.6136 | 0.5417 | 75.83 |
| 554 | 19.1200 | 0.5846 | 19.6994 | 0.6150 | 0.5434 | 75.92 |
| 100 | 19.4733 | 0.5518 | 19.9456 | 0.5823 | 0.5333 | 73.50 |
| 50 | 19.7954 | 0.5221 | 20.3229 | 0.5667 | 0.5045 | 72.67 |
| 30 | 20.2424 | 0.4819 | 20.9145 | 0.5417 | 0.4703 | 70.52 |
| 20 | 20.5885 | 0.4448 | 21.5472 | 0.4990 | 0.4272 | 68.75 |
| 10 | 21.8366 | 0.3384 | 22.6816 | 0.3700 | 0.3619 | 62.12 |

- Best CV rmse and TEST rmse are at n=644; best TEST exp_med is at n=554 (0.6150).
- TEST exp_med stays near 0.60 down to about n=200, then falls; it first drops below OligoAI (0.4349) at n=14.
- The 10 features left at n=10: `expr_target_dom_fraction`, `fold_mfe_aso5end`, `hybr_dna_dna_minus_dna_rna_dg`, `rbp_igf2bp2_aff_5`, `selfaso_homodimer_mismatched_dg`, `seq_gc_content`, `shape_slide_wing5`, `structure_sense_dist_to_splice_junction_exonic`, `structure_sense_signed_dist_to_canonical_start`, `structure_sense_start_from_end`.
