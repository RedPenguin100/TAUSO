"""SHAP attribution for the deployed booster: what drives the ranking, by family and by feature.

Panel a sums each ASO's SHAP values within a feature family and averages the magnitude over
ASOs, so a family that pulls in both directions on the same ASO reports the net pull rather
than the sum of its parts. Panel b keeps the features separate: one row per feature, one point
per ASO, placed at its SHAP value and coloured by where that ASO's feature value falls in the
distribution -- so a row reads as effect size, direction and the value that produced it.

Attribution is exact TreeSHAP, from the booster's own `pred_contribs`, over the held-out test
split. Writes the family table as a report, and the figure beside it.

    python notebooks/models/shap_figure.py [--version v1] [--top 16]
"""

import argparse
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import xgboost as xgb

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgba
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from notebooks.models import common

from tauso.inference.scoring import DEFAULT_VERSION, load_model

# Families follow the feature table published with the article; ASO:RNA shape is the duplex
# geometry block, which that table predates.
FAMILIES = (
    "Transcript context",
    "Sequence & motifs",
    "Target accessibility",
    "RBP binding",
    "Hybridization",
    "RNase H1 motifs",
    "Delivery & assay",
    "Translational context",
    "Off-target",
    "Expression & half-life",
    "Chemistry",
    "Toxicity / liability",
    "ASO:RNA shape",
)

# Hues evenly spaced in OKLCh at two lightness tiers, ordered so neighbours in this list stay
# apart under simulated colour-vision deficiency. Identity is carried by the labels; colour
# repeats it.
FAMILY_COLOUR = dict(
    zip(
        FAMILIES,
        (
            "#B45054", "#50BD86", "#7D62B7", "#D39838", "#0086A2", "#E281AD", "#598327",
            "#84A0F6", "#B05926", "#00BEB1", "#9A579E", "#B4A838", "#077CB9",
        ),
    )
)

VALUE_CMAP = LinearSegmentedColormap.from_list("feature_value", ["#2B6CB0", "#C9CDD4", "#C2352B"])
MISSING_COLOUR = "#D8D8D8"


def family_of(feature):
    """The family a feature column belongs to."""
    f = feature.lower()
    if f.startswith("rbp_"):
        return "RBP binding"
    if f.startswith("off_target"):
        return "Off-target"
    if f.startswith("rnase_"):
        return "RNase H1 motifs"
    if f.startswith("tox_"):
        return "Toxicity / liability"
    if f.startswith("shape_"):
        return "ASO:RNA shape"
    if f.startswith("hybr_") or f.startswith("on_target_total_hybridization"):
        return "Hybridization"
    if f.startswith(("structure_", "struct_", "flank_", "on_target_", "sense_")):
        return "Transcript context"
    if f.startswith(("fold_", "access_")):
        return "Target accessibility"
    if f.startswith(("cai_", "tai_", "enc_")):
        return "Translational context"
    if f.startswith(("chem_", "mod_")):
        return "Chemistry"
    if f.startswith(("expr_", "halflife")):
        return "Expression & half-life"
    if f.startswith(("transfection_", "density_", "volume_", "interaction_")):
        return "Delivery & assay"
    if f.startswith(("seq_", "ohe_", "dinuc_", "selfaso_")):
        return "Sequence & motifs"
    raise KeyError(f"No family for {feature!r}; add it to family_of().")


LABEL = {
    "interaction_selfaso_homodimer_perfect_dg_gymnosis": "ASO homodimer ΔG × gymnosis",
    "seq_internal_fold_rna": "ASO self-fold ΔG (RNA)",
    "sense_internal_fold": "Target self-fold ΔG",
    "selfaso_hp_dg": "ASO hairpin ΔG",
    "selfaso_homodimer_perfect_dg": "ASO homodimer ΔG",
    "seq_gc_content": "GC content (overall)",
    "seq_gc_content_3prime_end": "GC content (3′ end)",
    "mod_sugar_wing_gap_gc_delta": "Wing vs gap GC difference",
    "mod_sugar_gap_gc_content": "GC content of the gap",
    "hybr_cet_wing5_dg": "cEt 5′-wing duplex ΔG",
    "hybr_cet_wing3_dg": "cEt 3′-wing duplex ΔG",
    "hybr_dna_rna_dg": "On-target duplex ΔG (DNA:RNA)",
    "hybr_moe_md_gb_dg": "MOE duplex ΔG (MD, GB)",
    "expr_rnase_transcript": "RNase H1 expression",
    "expr_target": "Target mRNA expression",
    "expr_target_dom_fraction": "Canonical isoform share",
    "halflife_value": "Target mRNA half-life",
    "structure_sense_dist_to_splice_junction_exonic": "Distance to splice junction (exonic)",
    "structure_sense_dist_to_closest_splice_junction": "Distance to nearest splice junction",
    "structure_sense_start_norm": "Site position along the transcript",
    "transfection_gymnosis": "Delivery: gymnosis",
    "transfection_lipofection": "Delivery: lipofection",
    "transfection_electroporation": "Delivery: electroporation",
    "rnase_score_dinucleotide_R4a_dinuc_dynamic": "RNase-H R4a dinucleotide score",
    "rnase_krel_dinucleotide_score_R4a_krel_dinuc_dynamic": "RNase-H R4a dinucleotide k_rel score",
    "fold_mfe_win25_flank30_step4": "Target MFE (win25·flank30·step4)",
    "fold_mfe_win40_flank30_step5": "Target MFE (win40·flank30·step5)",
    "ohe_pos0_G": "5′-terminal base = G",
    "ohe_pos0_A": "5′-terminal base = A",
    "structure_sense_junction_signed_logdist_intronic": "Splice-junction distance, intronic (signed log)",
    "structure_sense_junction_logdist_exonic": "Splice-junction distance, exonic (log)",
    "structure_sense_signed_dist_to_canonical_start": "Distance to the start codon (signed)",
    "structure_sense_branch_point_signed_dist": "Offset from the branch point (nt)",
    "flank_gc_content_20": "Target flank GC content (±20 nt)",
    "seq_entropy": "Sequence entropy",
    "selfaso_hp_dg_per_nt": "ASO hairpin ΔG per nucleotide",
    "fold_mfe_aso5end": "Target MFE at the ASO 5′ end",
    "fold_mfe_aso3end": "Target MFE at the ASO 3′ end",
}


def label_of(feature):
    return LABEL.get(feature, feature)


def shap_matrix(booster, frame, features, cache):
    """Exact TreeSHAP for `frame`, one column per feature, cached on disk."""
    if cache.exists():
        stored = np.load(cache, allow_pickle=True)
        if list(stored["features"]) == list(features) and stored["values"].shape[0] == len(frame):
            return stored["values"]
    matrix = xgb.DMatrix(frame[features].to_numpy(np.float64), feature_names=features)
    values = booster.predict(matrix, pred_contribs=True)[:, :-1]  # the last column is the bias
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, values=values, features=np.array(features))
    return values


def family_importance(values, features):
    """Mean |sum of SHAP within the family| per ASO, and the family's feature count."""
    rows = []
    for family in FAMILIES:
        columns = [i for i, f in enumerate(features) if family_of(f) == family]
        if not columns:
            continue
        rows.append(
            {
                "family": family,
                "importance": float(np.abs(values[:, columns].sum(axis=1)).mean()),
                "n_features": len(columns),
            }
        )
    return pd.DataFrame(rows).sort_values("importance", ascending=False).reset_index(drop=True)


def _swarm_offsets(x, rng, height=0.42, bins=96):
    """Vertical offsets that spread points by local density, so a row reads as a distribution."""
    offsets = np.zeros(len(x))
    if len(x) == 0 or not np.isfinite(x).any() or np.nanmax(x) == np.nanmin(x):
        return offsets
    edges = np.linspace(np.nanmin(x), np.nanmax(x), bins + 1)
    which = np.clip(np.digitize(x, edges) - 1, 0, bins - 1)
    counts = np.bincount(which, minlength=bins)
    busiest = max(counts.max(), 1)
    for b in range(bins):
        members = np.flatnonzero(which == b)
        if len(members) == 0:
            continue
        spread = height * np.sqrt(len(members) / busiest)
        step = np.linspace(-spread, spread, len(members)) if len(members) > 1 else np.zeros(1)
        offsets[members] = rng.permutation(step)
    return offsets


def _colour_by_value(series):
    """Point colours from the feature's own distribution, by percentile so outliers can't flatten it."""
    values = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
    known = np.isfinite(values)
    shade = np.full(len(values), np.nan)
    if known.any():
        shade[known] = pd.Series(values[known]).rank(pct=True).to_numpy()
    colours = np.tile(to_rgba(MISSING_COLOUR), (len(values), 1))
    colours[known] = VALUE_CMAP(shade[known])
    return colours


def draw(families, values, features, frame, top, out_path, subtitle):
    order = np.argsort(-np.abs(values).mean(axis=0))[:top]
    rng = np.random.default_rng(0)

    figure = plt.figure(figsize=(10.5, 4.2 + 0.46 * top))
    grid = figure.add_gridspec(
        3, 1, height_ratios=[len(families) * 0.30, top * 0.46, 1.0], hspace=0.38,
    )

    ax_family = figure.add_subplot(grid[0])
    bars = ax_family.barh(
        families["family"], families["importance"],
        color=[FAMILY_COLOUR[f] for f in families["family"]], height=0.72,
    )
    ax_family.invert_yaxis()
    ax_family.set_xlim(0, families["importance"].max() * 1.26)
    for bar, importance, count in zip(bars, families["importance"], families["n_features"]):
        ax_family.text(
            bar.get_width() + families["importance"].max() * 0.015, bar.get_y() + bar.get_height() / 2,
            f"{importance:.2f}  ({count} feat.)", va="center", ha="left", fontsize=8, color="#333333",
        )
    ax_family.set_xlabel(
        "family importance = mean | Σ SHAP | (percentage points of inhibition)\n"
        "Σ summed over the family's features per ASO — within-family effects can cancel",
        fontsize=8.5,
    )
    ax_family.set_title("a   What drives the ranking — by feature family", loc="left", fontsize=11, weight="bold")
    ax_family.tick_params(axis="y", length=0, labelsize=9)
    for side in ("top", "right", "left"):
        ax_family.spines[side].set_visible(False)
    ax_family.grid(axis="x", color="#E6E6E6", lw=0.6)
    ax_family.set_axisbelow(True)

    ax_beeswarm = figure.add_subplot(grid[1])
    for row, column in enumerate(order):
        feature = features[column]
        shap_column = values[:, column]
        y = (top - 1 - row) + _swarm_offsets(shap_column, rng)
        ax_beeswarm.scatter(
            shap_column, y, s=2.0, c=_colour_by_value(frame[feature]), linewidths=0, alpha=0.65, rasterized=True,
        )
        # Between the label and the plotting area, so it never collides with either.
        ax_beeswarm.add_patch(
            plt.Rectangle(
                (-0.020, top - 1 - row - 0.28), 0.012, 0.56,
                transform=ax_beeswarm.get_yaxis_transform(), clip_on=False,
                color=FAMILY_COLOUR[family_of(feature)],
            )
        )

    ax_beeswarm.axvline(0, color="#9AA0A6", lw=0.9)
    ax_beeswarm.set_ylim(-0.7, top - 0.3)
    ax_beeswarm.set_yticks(range(top))
    ax_beeswarm.set_yticklabels(
        [f"#{rank + 1}  {label_of(features[c])}" for rank, c in reversed(list(enumerate(order)))], fontsize=9,
    )
    ax_beeswarm.tick_params(axis="y", length=0, pad=14)
    ax_beeswarm.set_xlabel("SHAP value (impact on predicted inhibition, percentage points)", fontsize=9)
    ax_beeswarm.set_title("b   Top features — effect, direction & value", loc="left", fontsize=11, weight="bold")
    for side in ("top", "right", "left"):
        ax_beeswarm.spines[side].set_visible(False)
    ax_beeswarm.grid(axis="x", color="#EFEFEF", lw=0.6)
    ax_beeswarm.set_axisbelow(True)

    bar = figure.colorbar(
        plt.cm.ScalarMappable(norm=Normalize(0, 1), cmap=VALUE_CMAP), ax=ax_beeswarm, pad=0.01, fraction=0.025,
    )
    bar.set_ticks([0, 1])
    bar.set_ticklabels(["low", "high"])
    bar.set_label("feature value (percentile)", fontsize=8.5)
    bar.outline.set_visible(False)

    ax_legend = figure.add_subplot(grid[2])
    ax_legend.axis("off")
    ax_legend.legend(
        handles=[Patch(facecolor=FAMILY_COLOUR[f], label=f) for f in families["family"]],
        loc="upper center", ncol=4, frameon=False, fontsize=8.5, title="feature family (the swatches in b)",
        title_fontsize=8.5,
    )

    figure.suptitle(
        f"TAUSO feature attribution — {subtitle}", x=0.02, ha="left", fontsize=12.5, weight="bold",
    )
    figure.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--version", default=DEFAULT_VERSION, help="deployed model version to attribute")
    ap.add_argument("--top", type=int, default=16, help="features shown in panel b")
    args = ap.parse_args()

    booster, features = load_model(args.version)
    df, _ = common.load_dataset()
    _, test = common.split(df)

    values = shap_matrix(booster, test, features, common.RESULTS_DIR / f"shap_{args.version}_test.npz")
    families = family_importance(values, features)

    lines = [
        f"Deployed model {args.version}: SHAP over the held-out test split, n={len(test)}",
        f"{len(features)} features in {len(families)} families",
        "",
        f"{'family':24}{'mean |sum SHAP|':>17}{'features':>10}",
        *(f"{r.family:24}{r.importance:17.3f}{r.n_features:10d}" for r in families.itertuples()),
    ]
    common.write_report(lines, f"shap_{args.version}_family_importance.txt")

    out_path = common.RESULTS_DIR / f"shap_{args.version}.png"
    draw(families, values, features, test, args.top, out_path, f"deployed model {args.version}, held-out test split")
    print(f"saved -> {out_path}")


if __name__ == "__main__":
    main()
