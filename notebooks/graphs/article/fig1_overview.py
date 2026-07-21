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
DELIV_ZOOM = 0.085        # thumbnail size of the delivery insets
OFFT_IMG = ASSETS / "offtarget.png"          # off-target: target vs off-target gene + "?"
OFFT_ZOOM = 0.10          # thumbnail size of the off-target inset
ASO_IMG = ASSETS / "seq_aso.png"             # lettered gapmer strand for sequence composition & motifs
ASO_ZOOM = 0.17           # thumbnail size of the sequence strand
SEQ_CLOUDS = ["Sequence composition\n(GC · AT-skew)", "RNase-H1 motif", "RBP motif"]

TAUSO_CALL = "TAUSO learns these signals\nand ranks candidate ASOs\nby predicted knockdown"

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


def leader(ax, xy1, xy2, color=GREY, lw=1.0):
    ax.plot([xy1[0], xy2[0]], [xy1[1], xy2[1]], color=color, lw=lw, zorder=3, solid_capstyle="round")


def inset_img(ax, path, xy, zoom, z=6, align=(0.5, 0.0)):
    ax.add_artist(AnnotationBbox(OffsetImage(plt.imread(str(path)), zoom=zoom), xy,
                                 frameon=False, box_alignment=align, zorder=z))


def cloud(ax, x, y, text, fs=8.4):
    ax.text(x, y, text, fontsize=fs, ha="center", va="center", color=INK, zorder=9,
            bbox=dict(boxstyle="round,pad=0.5", fc="#eef2f7", ec="#9fb0c3", lw=1.1))


# ------------------------------------------------------------------ panel a
def panel_a(ax):
    ax.set_xlim(-4, 104); ax.set_ylim(5, 99); ax.axis("off")
    ax.text(-3, 97.5, "a", fontsize=19, fontweight="bold", va="top")
    ax.text(4, 97.5, "What determines gapmer ASO efficacy", fontsize=13.5, fontweight="bold", va="top", color=INK)

    BF, BE = "#f7f8fa", "#c9d2da"                              # showcase box fill / edge

    # ===== TOP ROW: Hybridization | Off-target | Accessibility =====
    rbox(ax, 1, 77, 29, 17.5, fc=BF, ec=BE, lw=1.1, rad=1.6, z=0)
    ax.text(15.5, 92.4, "Hybridization affinity", fontsize=9.6, ha="center", va="center", fontweight="bold", color=INK)
    inset_img(ax, HYBR_IMGS[0], (15.5, 87), HYBR_ZOOM, align=(0.5, 0.5))
    inset_img(ax, HYBR_IMGS[1], (15.5, 81), HYBR_ZOOM, align=(0.5, 0.5))

    rbox(ax, 33.5, 77, 31, 17.5, fc=BF, ec=BE, lw=1.1, rad=1.6, z=0)
    ax.text(49, 92.4, "Off-target", fontsize=9.6, ha="center", va="center", fontweight="bold", color=INK)
    inset_img(ax, OFFT_IMG, (49, 83.5), OFFT_ZOOM, align=(0.5, 0.5))
    ax.text(36.5, 90.2, "✗", fontsize=15, ha="center", va="center", color="#C0392B", fontweight="bold", zorder=8)
    ax.text(61.5, 90.2, "✓", fontsize=15, ha="center", va="center", color=GREEN, fontweight="bold", zorder=8)

    rbox(ax, 68, 77, 30, 17.5, fc=BF, ec=BE, lw=1.1, rad=1.6, z=0)
    ax.text(83, 92.4, "Accessibility", fontsize=9.6, ha="center", va="center", fontweight="bold", color=INK)
    inset_img(ax, ACCESS_IMGS[0], (77, 79), ACCESS_ZOOM)
    inset_img(ax, ACCESS_IMGS[1], (89, 79), ACCESS_ZOOM)

    # ===== CENTRE: mechanism (cell -> ASO:target duplex -> RNase H1) + hand-off to TAUSO =====
    rbox(ax, 8, 53, 52, 20, fc=BLUE, ec="#bcd3e2", lw=1.2, rad=3.0, alpha=0.10, z=0)
    yT = 63.0
    for (x0, x1) in [(12, 22), (27, 43), (48, 58)]:
        rbox(ax, x0, yT - 1.3, x1 - x0, 2.6, fc=BLUE, ec=BLUE, rad=0.6, z=2)
    for (a, b) in [(22, 27), (43, 48)]:
        ax.plot([a, b], [yT, yT], color=GREY, lw=1.4, zorder=1)
    ax.text(12, yT - 3.1, "pre-mRNA target", fontsize=8.2, ha="left", va="center", color="#5b6b78")
    aso_x0, aso_x1, yA = 28, 42, 67.6
    rbox(ax, aso_x0, yA - 1.2, aso_x1 - aso_x0, 2.4, fc=ACCENT, ec=ACCENT, rad=0.6, z=3)
    for i in range(7):
        xt = aso_x0 + 1.6 + i * 1.8
        ax.plot([xt, xt], [yT + 1.3, yA - 1.2], color="#6b7885", lw=1.3, zorder=2)   # base-pairing (between the RNAs)
    ax.text(35, yA + 1.8, "gapmer ASO", fontsize=9.6, ha="center", va="bottom", fontweight="bold", color=ACCENT)
    ax.add_patch(Wedge((35, 59.0), 3.5, 122, 58, fc=PURPLE, ec=PURPLE, lw=1.0, alpha=0.85, zorder=2))  # RNase H1 (pac-man)
    ax.text(35, 53.9, "RNase H1", fontsize=8.2, ha="center", va="center", color=PURPLE, fontweight="bold")
    arrow(ax, (35, 60.9), (35, 61.9), color=PURPLE, lw=1.6, ms=11)          # cleavage
    arrow(ax, (60, 63), (66, 63), color=INK, lw=2.4, ms=18)
    rbox(ax, 67, 56, 33, 14, fc=ACCENT_L, ec=ACCENT, lw=1.8, rad=1.6, z=2)
    ax.text(83.5, 63, TAUSO_CALL, fontsize=9.4, ha="center", va="center", color=INK, zorder=6)

    # ===== DELIVERY ROW =====
    rbox(ax, 14, 32, 72, 18, fc=BF, ec=BE, lw=1.1, rad=1.6, z=0)
    ax.text(50, 47.6, "Delivery / uptake", fontsize=9.6, ha="center", va="center", fontweight="bold", color=INK)
    for img, cx in zip(DELIV_IMGS, (33, 50, 67)):
        inset_img(ax, img, (cx, 40), DELIV_ZOOM, align=(0.5, 0.5))

    # ===== SEQUENCE COMPOSITION & MOTIFS (full width): lettered strand + extracted 'clouds' =====
    rbox(ax, 1, 6, 97, 23, fc=BF, ec=BE, lw=1.1, rad=1.6, z=0)
    ax.text(49.5, 26.6, "Sequence composition & motifs", fontsize=9.6, ha="center", va="center", fontweight="bold", color=INK)
    inset_img(ax, ASO_IMG, (49.5, 11.5), ASO_ZOOM, align=(0.5, 0.5))
    for (cx, anchor) in [(20, 27), (49.5, 48), (79, 71)]:
        leader(ax, (cx, 18.3), (anchor, 13.6), color="#b9c2cb", lw=0.9)
        ax.scatter([anchor], [13.6], s=8, color="#8a94a0", zorder=7)
    cloud(ax, 20, 21, SEQ_CLOUDS[0])
    cloud(ax, 49.5, 21.2, SEQ_CLOUDS[1])
    cloud(ax, 79, 21, SEQ_CLOUDS[2])


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
    rbox(ax, 5, 10, 24, 40, fc="white", ec=BLUE, lw=1.1, rad=1.2, z=2)
    ax.text(cx, 47, "Features per ASO", fontsize=9.4, ha="center", va="top", fontweight="bold", color=BLUE)
    for i, f in enumerate(FEATURE_CLASSES):
        ax.text(7.5, 42.5 - i * 5.4, "•", fontsize=10, ha="left", va="top", color=BLUE)
        ax.text(9.5, 42.5 - i * 5.4, f, fontsize=8.7, ha="left", va="top", color=INK)

    # ---- MODEL column
    mx = 50
    for i, line in enumerate(MODEL_LINES):
        y = 74 - i * 12.5
        strong = i == 0
        rbox(ax, 38, y - 4.4, 24, 8.8, fc="white", ec=GREEN if strong else "#bcd6c2",
             lw=1.5 if strong else 1.0, rad=1.2, z=2)
        ax.text(mx, y, line, fontsize=9.2 if strong else 8.7, ha="center", va="center",
                color=INK, fontweight="bold" if strong else "normal", zorder=3)
        if i < len(MODEL_LINES) - 1:
            ax.plot([mx, mx], [y - 4.4, y - 8.1], color="#bcd6c2", lw=1.2, zorder=1)

    # ---- APPLICATION column
    ax_ = 83
    for i, line in enumerate(APP_LINES):
        y = 74 - i * 12
        rbox(ax, 71, y - 4.4, 24, 8.8, fc="white", ec=AMBER, lw=1.2, rad=1.2, z=2)
        ax.text(ax_, y, line, fontsize=8.9, ha="center", va="center", color=INK, zorder=3)
    # MD as a distinct, optional (dashed) downstream refinement stage
    rbox(ax, 71, 22, 24, 9.5, fc="white", ec=PURPLE, lw=1.4, rad=1.2, ls=(0, (4, 2)), z=2)
    ax.text(ax_, 26.8, MD_STAGE, fontsize=8.9, ha="center", va="center", color=PURPLE, fontweight="bold", zorder=3)
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
