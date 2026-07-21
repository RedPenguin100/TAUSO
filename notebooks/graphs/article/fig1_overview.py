"""Fig 1 — study overview (DRAFT).

a) Conceptual schematic: a gapmer ASO hybridizes its target transcript and directs RNase H1 cleavage; the
   determinants of knockdown efficacy (hybridization, accessibility, off-target, RNase-H1, delivery, toxicity)
   are the signals TAUSO learns to rank candidate ASOs.
b) Pipeline: Data (ASO Atlas -> filtered gapmers -> per-ASO feature classes) -> Model (XGBoost, within-
   experiment ranking objective, gene-grouped CV, per-ASO features) -> Application (rank / design / benchmark;
   MD refinement of top candidates).

Panel-a artwork lives beside this file in fig1_assets/ (self-contained; no external paths).

Run: python fig1_overview.py
"""
import sys
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Ellipse, Arc, Circle, Wedge
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # graphs/ for consts
from consts import save, style, ACCENT, BLUE, GREY, INK

# ================================= TWEAK ME =================================
SOURCE   = "ASO Atlas"          # data source name (matches the text; NOT "OligoAtlas")
RAW_N    = "178,624"            # raw ASO–target knockdown measurements
FILT_N   = "143,904"            # after filtering to strict gapmers (v15 averaged table)
N_FEAT   = "485"                # deployed per-ASO features (all non-zero-variance)
GENE_OBJ = "within-experiment ranking"   # the clean_exp objective, in plain words

# panel a — accessibility showcase: open (accessible) vs closed (inaccessible) target illustrations
ASSETS = Path(__file__).resolve().parent / "fig1_assets"
ACCESS_IMGS = (ASSETS / "accessible_target.png", ASSETS / "inaccessible_target.png")
ACCESS_ZOOM = 0.072       # thumbnail size of the accessibility insets
HYBR_IMGS = (ASSETS / "weak_hybridization.png", ASSETS / "hybr_strong.png")   # weak vs strong duplex
HYBR_ZOOM = 0.128         # thumbnail size of the hybridization insets
DELIV_IMGS = (ASSETS / "delivery_lipofection.png", ASSETS / "delivery_electroporation.png",
              ASSETS / "delivery_freeuptake.png")   # lipofection · electroporation · free uptake
DELIV_ZOOM = 0.098        # thumbnail size of the delivery insets
DELIV_LABELS = ["lipofection", "electroporation", "gymnosis"]   # gymnosis = free (gymnotic) uptake
OFFT_IMG = ASSETS / "offtarget.png"          # off-target: target vs off-target gene + "?"
OFFT_ZOOM = 0.10          # thumbnail size of the off-target inset
ASO_IMG = ASSETS / "seq_aso.png"             # lettered gapmer strand for sequence composition & motifs
ASO_ZOOM = 0.166          # thumbnail size of the sequence strand (spans the row)

# panel b — feature classes (per-ASO), model card, applications
FEATURE_CLASSES = [
    "Sequence & chemistry", "Hybridization thermodynamics", "Target accessibility",
    "Off-target burden", "RBP / protein context", "Toxicity / liability",
]
MODEL_LINES = [
    "XGBoost gradient-boosted trees",
    f"Objective: {GENE_OBJ}",
    "Gene-grouped cross-validation",
    f"{N_FEAT} per-ASO features",
    "Tuned, regularized parameters",
]
APP_LINES = [
    "Rank & prioritize candidate ASOs",
    "Design new ASOs (e.g. DIAPH3)",
    "Benchmark vs OligoAI",
]
MD_STAGE = "MD refinement of top candidates"   # shown as a distinct downstream stage (dashed = optional)

# colours
GREEN, AMBER, PURPLE = "#3C8C4E", "#C89A2B", "#7B5EA7"
BLUE_L, GREEN_L, AMBER_L, ACCENT_L, GREY_L = "#E7F1F7", "#E9F3EB", "#FAF1DF", "#FBEAE3", "#EEF1F3"
# ===========================================================================


def rbox(ax, x, y, w, h, fc="none", ec=INK, lw=1.3, rad=1.4, ls="-", alpha=1.0, z=1):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={rad}",
                                fc=fc, ec=ec, lw=lw, ls=ls, alpha=alpha, zorder=z, mutation_aspect=1.0))


def arrow(ax, xy1, xy2, color=INK, lw=2.2, ms=16, ls="-", z=5):
    ax.add_patch(FancyArrowPatch(xy1, xy2, arrowstyle="-|>", mutation_scale=ms, lw=lw,
                                 color=color, ls=ls, shrinkA=0, shrinkB=0, zorder=z))


def inset_img(ax, path, xy, zoom, z=6, align=(0.5, 0.0)):
    ax.add_artist(AnnotationBbox(OffsetImage(plt.imread(str(path)), zoom=zoom), xy,
                                 frameon=False, box_alignment=align, zorder=z))


# ------------------------------------------------------------------ panel a
def panel_a(ax):
    ax.set_xlim(-4, 104); ax.set_ylim(1, 99); ax.axis("off")
    ax.text(-3, 97.5, "a", fontsize=19, fontweight="bold", va="top")
    ax.text(4, 97.5, "What determines gapmer ASO efficacy", fontsize=13.5, fontweight="bold", va="top", color=INK)

    BF, BE = "#f7f8fa", "#c9d2da"                              # showcase box fill / edge

    # ===== ROW 1: mechanism (cell -> ASO:target duplex -> RNase H1) =====
    rbox(ax, 24, 73, 52, 20, fc=BLUE, ec="#bcd3e2", lw=1.2, rad=3.0, alpha=0.10, z=0)
    yT = 83.0
    for (x0, x1) in [(28, 38), (43, 59), (64, 74)]:
        rbox(ax, x0, yT - 1.3, x1 - x0, 2.6, fc=BLUE, ec=BLUE, rad=0.6, z=2)
    for (a, b) in [(38, 43), (59, 64)]:
        ax.plot([a, b], [yT, yT], color=GREY, lw=1.4, zorder=1)
    ax.text(28, yT - 3.1, "pre-mRNA target", fontsize=8.2, ha="left", va="center", color="#5b6b78")
    aso_x0, aso_x1, yA = 44, 58, 87.6
    rbox(ax, aso_x0, yA - 1.2, aso_x1 - aso_x0, 2.4, fc=ACCENT, ec=ACCENT, rad=0.6, z=3)
    for i in range(7):
        xt = aso_x0 + 1.6 + i * 1.8
        ax.plot([xt, xt], [yT + 1.3, yA - 1.2], color="#6b7885", lw=1.3, zorder=2)   # base-pairing (between the RNAs)
    ax.text(51, yA + 1.8, "gapmer ASO", fontsize=9.6, ha="center", va="bottom", fontweight="bold", color=ACCENT)
    ax.add_patch(Wedge((51, 79.0), 3.5, 122, 58, fc=PURPLE, ec=PURPLE, lw=1.0, alpha=0.85, zorder=2))  # RNase H1 (pac-man)
    ax.text(51, 74.0, "RNase H1", fontsize=8.2, ha="center", va="center", color=PURPLE, fontweight="bold")
    arrow(ax, (51, 80.9), (51, 81.9), color=PURPLE, lw=1.6, ms=11)          # cleavage

    # ===== ROW 2: determinants — Hybridization | Off-target | Accessibility =====
    rbox(ax, 1, 51, 29, 17.5, fc=BF, ec=BE, lw=1.1, rad=1.6, z=0)
    ax.text(15.5, 66.4, "Hybridization affinity", fontsize=9.6, ha="center", va="center", fontweight="bold", color=INK)
    inset_img(ax, HYBR_IMGS[0], (15.5, 61), HYBR_ZOOM, align=(0.5, 0.5))
    inset_img(ax, HYBR_IMGS[1], (15.5, 55), HYBR_ZOOM, align=(0.5, 0.5))

    rbox(ax, 33.5, 51, 31, 17.5, fc=BF, ec=BE, lw=1.1, rad=1.6, z=0)
    ax.text(49, 66.4, "Off-target", fontsize=9.6, ha="center", va="center", fontweight="bold", color=INK)
    inset_img(ax, OFFT_IMG, (49, 57.5), OFFT_ZOOM, align=(0.5, 0.5))
    ax.text(36.5, 64.2, "✗", fontsize=15, ha="center", va="center", color="#C0392B", fontweight="bold", zorder=8)
    ax.text(61.5, 64.2, "✓", fontsize=15, ha="center", va="center", color=GREEN, fontweight="bold", zorder=8)

    rbox(ax, 68, 51, 30, 17.5, fc=BF, ec=BE, lw=1.1, rad=1.6, z=0)
    ax.text(83, 66.4, "Accessibility", fontsize=9.6, ha="center", va="center", fontweight="bold", color=INK)
    inset_img(ax, ACCESS_IMGS[0], (77, 53), ACCESS_ZOOM)
    inset_img(ax, ACCESS_IMGS[1], (89, 53), ACCESS_ZOOM)

    # ===== ROW 3: delivery — one box per method (same x-extent as the determinant row) =====
    for (x0, x1), img, lbl in zip([(1, 30), (33.5, 64.5), (68, 98)], DELIV_IMGS, DELIV_LABELS):
        cx = (x0 + x1) / 2
        rbox(ax, x0, 27, x1 - x0, 20, fc=BF, ec=BE, lw=1.1, rad=1.6, z=0)
        ax.text(cx, 44.5, lbl.capitalize(), fontsize=9.6, ha="center", va="center", fontweight="bold", color=INK)
        inset_img(ax, img, (cx, 35.0), DELIV_ZOOM, align=(0.5, 0.5))

    # ===== ROW 4: sequence composition & motifs — full-width strand -> feature calculations =====
    rbox(ax, 1, 3, 97, 21, fc=BF, ec=BE, lw=1.1, rad=1.6, z=0)
    ax.text(49.5, 21.6, "Sequence composition & motifs", fontsize=9.6, ha="center", va="center", fontweight="bold", color=INK)
    inset_img(ax, ASO_IMG, (49, 8.5), ASO_ZOOM, align=(0.5, 0.5))              # strand spans the row
    rbox(ax, 66, 14.6, 30, 5.0, fc="white", ec=BLUE, lw=1.2, rad=1.3, z=4)     # feature-calculations node (upper-right)
    ax.text(81, 17.1, "Feature calculations", fontsize=9.0, ha="center", va="center", fontweight="bold", color=BLUE, zorder=5)
    arrow(ax, (81, 12.7), (81, 14.4), color=GREY, lw=1.6, ms=12)              # strand -> feature calculations


# ------------------------------------------------------------------ panel b
def panel_b(ax):
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
    ax.text(1, 99, "b", fontsize=19, fontweight="bold", va="top")
    ax.text(6, 98, "The TAUSO pipeline", fontsize=13.5, fontweight="bold", va="top", color=INK)

    cols = [(3, 31, "DATA", BLUE, BLUE_L), (36, 64, "MODEL", GREEN, GREEN_L),
            (69, 97, "APPLICATION", AMBER, AMBER_L)]
    for (x0, x1, head, c, cl) in cols:
        rbox(ax, x0, 5, x1 - x0, 82, fc=cl, ec=c, lw=1.5, rad=1.6, z=1)          # container
        rbox(ax, x0, 88, x1 - x0, 6.5, fc=c, ec=c, rad=1.2, z=2)                 # header bar
        ax.text((x0 + x1) / 2, 91.3, head, fontsize=12, ha="center", va="center",
                color="white", fontweight="bold", zorder=3)

    # inter-column flow arrows
    arrow(ax, (31, 46), (36, 46), color=INK, lw=2.4, ms=17)
    arrow(ax, (64, 46), (69, 46), color=INK, lw=2.4, ms=17)

    # ---- DATA column
    cx = 17
    ax.text(cx, 83, SOURCE, fontsize=11.5, ha="center", fontweight="bold", color=BLUE)
    ax.text(cx, 79.5, f"{RAW_N} ASO–target\nknockdown measurements", fontsize=8.8, ha="center", va="top", color=INK)
    arrow(ax, (cx, 73.5), (cx, 71), color=GREY, lw=1.6, ms=12)
    ax.text(cx, 70, "filter → strict gapmers", fontsize=9, ha="center", va="top", color="#5b6b78", style="italic")
    ax.text(cx, 66.5, f"{FILT_N} ASOs", fontsize=10.5, ha="center", va="top", fontweight="bold", color=INK)
    arrow(ax, (cx, 62.5), (cx, 60), color=GREY, lw=1.6, ms=12)
    ax.text(cx, 59, "gene-grouped\ntrain / test split", fontsize=8.8, ha="center", va="top", color=INK)
    rbox(ax, 4, 10, 26, 40, fc="white", ec=BLUE, lw=1.1, rad=1.2, z=2)
    ax.text(cx, 47, "Features per ASO", fontsize=9.4, ha="center", va="top", fontweight="bold", color=BLUE)
    for i, f in enumerate(FEATURE_CLASSES):
        ax.text(6.0, 42.5 - i * 5.4, "•", fontsize=10, ha="left", va="top", color=BLUE)
        ax.text(8.0, 42.5 - i * 5.4, f, fontsize=8.0, ha="left", va="top", color=INK)

    # ---- MODEL column
    mx = 50
    for i, line in enumerate(MODEL_LINES):
        y = 74 - i * 12.5
        strong = i == 0
        rbox(ax, 37, y - 4.4, 26, 8.8, fc="white", ec=GREEN if strong else "#bcd6c2",
             lw=1.5 if strong else 1.0, rad=1.2, z=2)
        ax.text(mx, y, line, fontsize=8.6 if strong else 8.2, ha="center", va="center",
                color=INK, fontweight="bold" if strong else "normal", zorder=3)
        if i < len(MODEL_LINES) - 1:
            ax.plot([mx, mx], [y - 4.4, y - 8.1], color="#bcd6c2", lw=1.2, zorder=1)

    # ---- APPLICATION column
    ax_ = 83
    for i, line in enumerate(APP_LINES):
        y = 74 - i * 12
        rbox(ax, 70, y - 4.4, 26, 8.8, fc="white", ec=AMBER, lw=1.2, rad=1.2, z=2)
        ax.text(ax_, y, line, fontsize=8.3, ha="center", va="center", color=INK, zorder=3)
    # MD as a distinct, optional (dashed) downstream refinement stage
    rbox(ax, 70, 22, 26, 9.5, fc="white", ec=PURPLE, lw=1.4, rad=1.2, ls=(0, (4, 2)), z=2)
    ax.text(ax_, 26.8, MD_STAGE, fontsize=8.2, ha="center", va="center", color=PURPLE, fontweight="bold", zorder=3)
    arrow(ax, (ax_, 33.6), (ax_, 31.5), color=PURPLE, lw=1.4, ms=11, ls="-")
    ax.text(ax_, 15, "prioritized, validated\nASO candidates", fontsize=8.6, ha="center", va="center",
            color="#7a5c12", style="italic", zorder=3)


style()
fig = plt.figure(figsize=(11.5, 16.5))
gs = fig.add_gridspec(2, 1, height_ratios=[1.6, 1.0], hspace=0.05)
panel_a(fig.add_subplot(gs[0]))
panel_b(fig.add_subplot(gs[1]))

CAPTION = f"""
**Figure 1. Study overview.**
(a) A gapmer antisense oligonucleotide (ASO) hybridizes to its target transcript and recruits RNase H1, which
cleaves the RNA:DNA duplex to knock down the target. Knockdown efficacy is governed by several competing
determinants — hybridization affinity, local target accessibility (RNA secondary structure), off-target
hybridization, RNase-H1 recruitment, cellular delivery/uptake, and sequence-related toxicity liabilities.
TAUSO learns these signals to rank candidate ASOs by predicted knockdown.
(b) Pipeline. Data: {RAW_N} ASO–target measurements from the {SOURCE}, filtered to {FILT_N} strict gapmers,
split by gene into train/test; each ASO is described by six feature classes. Model: XGBoost gradient-boosted
trees trained on a {GENE_OBJ} objective with gene-grouped cross-validation; {N_FEAT} per-ASO features and
tuned, regularized parameters define the final model. Application: ranking and design of candidate ASOs
(e.g. DIAPH3), benchmarking against OligoAI, and molecular-dynamics refinement of the top candidates.
"""
VEC = Path(__file__).resolve().parents[1] / "out" / "article"
VEC.mkdir(parents=True, exist_ok=True)
for ext in ("svg", "pdf"):                                    # editable vector copies (Inkscape / Illustrator)
    fig.savefig(VEC / f"fig1_overview.{ext}", bbox_inches="tight", transparent=True)
    print("wrote", (VEC / f"fig1_overview.{ext}").name)
save(fig, "article", "fig1_overview", caption=CAPTION)        # PNG (+ white copy + caption); closes fig
print("DONE")
