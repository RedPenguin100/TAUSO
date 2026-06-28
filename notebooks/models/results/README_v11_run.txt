TAUSO model-selection pipeline — v11 run (how to reproduce exactly)

STORE: oligo_features_v11.parquet | Zenodo 20981280 | md5 f056c919f5b791e501dccdc49d807f80
       (gymnosis interaction = internal_fold x gymnosis; 510 loader features)
CODE : branch feat/model-selection-pipeline @ commit 08453b3 (origin/main 7c1519e merged = v11)
ENV  : conda env tauso_research (Py3.11); tauso installed editable from TAUSO-min/src
DETERMINISM: frozen OligoAI split, GroupKFold (no RNG), fixed seeds (SEED/RANK_SEED=0, XGBoost seed=1),
             canonical column order -> run twice = identical reports.

Commands (run from /home/michael/career/tauso_article/TAUSO-min):
  python notebooks/models/1_parameter_choice.py                              # -> parameter_choice.txt
  python notebooks/models/2_objective_choice.py                             # -> objective_choice.txt
  python notebooks/models/3_selection_choice.py  --variant {clean_gene|clean_exp|pair_exp}
  python notebooks/models/4_feature_selection.py --variant {clean_gene|clean_exp|pair_exp}
