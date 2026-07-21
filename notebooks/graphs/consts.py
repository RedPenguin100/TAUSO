"""Shared output paths, palette, style and save helper for the graphs.

Code stays in graphs/<area>/...  ;  rendered figures go to graphs/out/<category>/
as <name>.png (transparent) plus a brief <name>.txt caption (paper figure legend).

Usage in a script:
    import sys; from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[N]))   # N = depth up to graphs/
    from consts import save, style, ACCENT, BLUE, GREY, INK
    style(); fig, ax = plt.subplots(...); ...
    save(fig, "raw_features/off_target", "off_target_burden_dose", caption="...")
"""
from pathlib import Path
import matplotlib.pyplot as plt
from PIL import Image

GRAPHS = Path(__file__).resolve().parent
OUT = GRAPHS / "out"


def _write_white(png_path):
    """Composite a transparent PNG onto a solid white background -> <name>_white.png."""
    im = Image.open(png_path).convert("RGBA")
    bg = Image.new("RGBA", im.size, (255, 255, 255, 255))
    white = png_path.with_name(png_path.stem + "_white.png")
    Image.alpha_composite(bg, im).convert("RGB").save(white)
    return white

ACCENT, BLUE, GREY, INK = "#E4572E", "#2E86AB", "#9AA5B1", "#1a1a1a"


def style():
    """House style with a transparent canvas (figure + axes)."""
    plt.rcParams.update({
        "font.family": "sans-serif", "font.size": 11,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": "#444",
        "figure.facecolor": "none", "axes.facecolor": "none", "savefig.facecolor": "none",
    })


def save(fig, category, name, caption=None, dpi=300):
    """Save a transparent PNG (and optional Markdown caption) under graphs/out/<category>/."""
    d = OUT / category
    d.mkdir(parents=True, exist_ok=True)
    fig.savefig(d / f"{name}.png", dpi=dpi, bbox_inches="tight", transparent=True)
    _write_white(d / f"{name}.png")                       # also a white-background copy
    if caption is not None:
        (d / f"{name}.md").write_text(caption.strip() + "\n")
    plt.close(fig)
    print("wrote", (d / f"{name}.png").relative_to(GRAPHS))
