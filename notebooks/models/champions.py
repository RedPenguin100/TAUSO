"""The six v15 HPO champions: the best trial per (objective x selection-metric) from the
distributed Optuna sweep on the v15 feature set (485 features, 121,997 train+val rows).

Exact, unrounded params. Each champion carries its own training objective (variant:
clean_exp / clean_gene / reg) and its tuned num_boost_round. cv_exp_med / cv_gxc_med are the
sweep's cross-validation scores for that trial, kept for provenance.

BACKEND holds only the execution knobs (tree method / device); set device="cpu" for CPU-only runs.
"""

BACKEND = dict(tree_method="hist", device="cuda")

# name -> variant (objective), tuned params, num_boost_round, and the sweep's CV scores.
CHAMPIONS = {
    "clean_exp_exp_med": {  # study clean_exp_v15, trial 229 (best CV exp_med)
        "variant": "clean_exp", "rounds": 6239, "cv_exp_med": 0.5503459783658935, "cv_gxc_med": 0.41274052071355527,
        "params": dict(max_depth=13, learning_rate=0.004998274438617547, min_child_weight=138,
                       subsample=0.7515264828254578, colsample_bytree=0.4182276770391009,
                       gamma=1.9369195974720954, reg_lambda=10.654538362400565, reg_alpha=2.9314301182576425)},
    "clean_exp_gxc_med": {  # study clean_exp_v15, trial 218 (best CV gxc_med)
        "variant": "clean_exp", "rounds": 5684, "cv_exp_med": 0.5476897787363348, "cv_gxc_med": 0.4275619858839182,
        "params": dict(max_depth=12, learning_rate=0.0070117697420548864, min_child_weight=128,
                       subsample=0.7934625082800734, colsample_bytree=0.40584567752603934,
                       gamma=4.428889791038342, reg_lambda=12.379192790873187, reg_alpha=2.49708380382469)},
    "clean_gene_exp_med": {  # study clean_gene_v15, trial 162 (best CV exp_med)
        "variant": "clean_gene", "rounds": 5941, "cv_exp_med": 0.5508666583162236, "cv_gxc_med": 0.4216414778364017,
        "params": dict(max_depth=12, learning_rate=0.004804277794768826, min_child_weight=131,
                       subsample=0.8490167887734463, colsample_bytree=0.5178869079192633,
                       gamma=5.550106202879606, reg_lambda=11.739195412555834, reg_alpha=0.7978672235057911)},
    "clean_gene_gxc_med": {  # study clean_gene_v15, trial 200 (best CV gxc_med)
        "variant": "clean_gene", "rounds": 6537, "cv_exp_med": 0.549281358457337, "cv_gxc_med": 0.4297783809840596,
        "params": dict(max_depth=14, learning_rate=0.004213552841413954, min_child_weight=130,
                       subsample=0.8488755408159513, colsample_bytree=0.41251148610000143,
                       gamma=5.243740580551425, reg_lambda=12.67076730219519, reg_alpha=1.7883859558155508)},
    "reg_exp_med": {  # study reg_v15, trial 11 (best CV exp_med)
        "variant": "reg", "rounds": 5662, "cv_exp_med": 0.540880014126899, "cv_gxc_med": 0.40264624758033535,
        "params": dict(max_depth=12, learning_rate=0.007304785337826039, min_child_weight=178,
                       subsample=0.9260494038493405, colsample_bytree=0.3052108502794771,
                       gamma=0.5964538891394238, reg_lambda=9.785168390391412, reg_alpha=0.9331332671603934)},
    "reg_gxc_med": {  # study reg_v15, trial 129 (best CV gxc_med)
        "variant": "reg", "rounds": 6902, "cv_exp_med": 0.516440860289358, "cv_gxc_med": 0.4083995283792484,
        "params": dict(max_depth=11, learning_rate=0.008702752339883183, min_child_weight=27,
                       subsample=0.8121111581278974, colsample_bytree=0.5870216308037386,
                       gamma=1.4310102483633413, reg_lambda=3.3923816009322296, reg_alpha=1.2496558456771425)},
}


def config(name):
    """Return (params, rounds, variant) for a champion, with the backend knobs merged into params."""
    c = CHAMPIONS[name]
    return {**BACKEND, **c["params"]}, c["rounds"], c["variant"]
