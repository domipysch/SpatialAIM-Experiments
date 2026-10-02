"""What every figure of the paper shares: the palette, the mapper names, and the
way a panel image is written.

Only what is genuinely common lives here. Each figure keeps its own rcParams and
its own geometry next to the code that uses them -- those are arguments about
that figure's layout, not house style.

The figures are assembled in Excalidraw: a panel is one PNG carrying only what is
anchored to data, and the letters, headings and legends are drawn around it. A
panel is therefore written at print size, and ``PanelWriter`` records that size
so the drawing can place the image at exactly the inches it was authored for.
"""

from __future__ import annotations

import logging
import struct
from pathlib import Path

import matplotlib as mpl
from matplotlib import pyplot as plt

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# geometry shared with the drawings
# ---------------------------------------------------------------------------
DPI = 500             # what a panel PNG is written at
SCALE = 300           # canvas units per final print inch, as in every panels.json
PREVIEW_DPI = 200     # the white-background copy kept for a quick look

# ---------------------------------------------------------------------------
# palette
# ---------------------------------------------------------------------------
INK = "#1e1e1e"
GREY = "#868e96"
PALE = "#DCDCDC"

# the three criteria, wherever a figure draws their scores
C_RECON = "#543005"
C_SPATIAL = "#C51B7D"
C_COHERENCE = "#3F007D"

# The cell-type families of figures 2.1, 2.3 and 2.4's sketches.
C_GLUT = "#D55E00"
C_GABA = "#0072B2"
C_NONN = "#009E73"

# The five spot-to-state mappers. All five hues are figure 2.1's; the two blues
# go to SpatialAIM's own parameter-free mappers, so that family is legible before the
# legend is consulted.
MAPPERS = ["nearest_centroid", "wann", "tangram", "tacco", "dot"]
MAPPER_COLOURS = {
    "nearest_centroid": "#0072B2",
    "wann": "#7fb8d8",
    "tangram": "#D55E00",
    "tacco": "#E69F00",
    "dot": "#009E73",
}
# How a mapper is written wherever a figure names it.
MAPPER_LABELS = {
    "nearest_centroid": "Nearest centroid",
    "wann": "WANN",
    "tangram": "Tangram",
    "tacco": "TACCO",
    "dot": "DOT",
}

# The agglomeration linkages, short because they share narrow columns.
TREES = ["ward", "average"]
TREE_LABELS = {"ward": "ward", "average": "avg"}

# Figure 3's heat maps: the plotting box inside each image as (left, bottom, width,
# height) fractions. Its panels.py draws the rows into it; its make_excalidraw.py
# places the linkage labels and mapper bars against those rows.
HEAT_AXES = (0.02, 0.265, 0.96, 0.475)

# The ten references in corpus order -- the column order of every heat map --
# and the short name each is labelled with.
REFERENCE_LABELS = {
    "01_mouse-ssp": "ssp",
    "02_colorectal-cancer": "colorectal",
    "03_mouse-liver": "liver",
    "04_mouse-kidney": "kidney",
    "05_mouse-embryo": "embryo",
    "06_mouse-retina": "retina",
    "07_mouse-mpfc": "mPFC",
    "08_mouse-visual-cortex": "visual",
    "09_human-breast-cancer": "breast",
    "10_human-mtg": "MTG",
}

# Difference in how many cell types a mapping used, on its own hue pair so that
# nobody reads it as a cosine: purple = SpatialAIM used fewer, green = more. Shared by
# figures 2.2 and 2.4, which both carry that quantity.
CMAP_DK = mpl.colors.LinearSegmentedColormap.from_list(
    "aim_dk", ["#3F007D", "#8C3CC8", "#C9A6E8", "#FFFFFF",
               "#A9DFCC", "#009E73", "#00543C"])
# Clipped, because MTG runs to +69 (76 annotated types against SpatialAIM's 9) and would
# flatten every other cell to white.
DK_LIMIT = 20


# The figures print at the paper's 5.17 in text width, so a 7.2 in drawing is
# scaled by 0.71 and every point size is authored 1/0.71 larger than it should
# print: 9.2 pt here comes out as 6.5 pt on the page. Figures 2.2 and 2.4 are
# authored at that size; figure 2.3's panel modules set their own, smaller base.
PAPER_RCPARAMS = {
    "font.size": 9.2, "axes.titlesize": 9.2, "axes.labelsize": 9.2,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "legend.fontsize": 8.0,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.0, "ytick.major.size": 2.0, "pdf.fonttype": 42,
}


def use_paper_rcparams(**overrides) -> None:
    """Apply ``PAPER_RCPARAMS``. Call it at import in a panel module: the sizes
    are part of how the panel is authored, so a figure that forgets it silently
    prints at matplotlib's defaults."""
    mpl.rcParams.update({**PAPER_RCPARAMS, **overrides})


def bare(ax) -> None:
    """Drop the two spines a scatter or box plot does not need."""
    ax.spines[["top", "right"]].set_visible(False)


def png_size(path: Path) -> tuple[int, int]:
    """Pixel size from the PNG header, to avoid a Pillow dependency."""
    with path.open("rb") as fh:
        head = fh.read(24)
    return struct.unpack(">II", head[16:24])


class PanelWriter:
    """Writes a figure's panels and remembers how large each one is.

    ``manifest`` is what ``make_excalidraw.py`` reads: per panel the file name
    and its size in canvas units, i.e. print inches times ``SCALE``. The drawing
    can then place every image at the size it was authored for instead of at
    whatever pixel count the DPI happened to produce.

    ``preview`` is an optional second copy on white at a screen resolution --
    handy to look at, never embedded anywhere.
    """

    def __init__(self, out: Path, preview: Path | None = None,
                 dpi: int = DPI, scale: int = SCALE):
        self.out = Path(out)
        self.preview = Path(preview) if preview is not None else None
        self.dpi = dpi
        self.scale = scale
        self.panels: dict[str, dict] = {}

    def save(self, fig, name: str) -> None:
        self.out.mkdir(parents=True, exist_ok=True)
        path = self.out / f"{name}.png"
        fig.savefig(path, dpi=self.dpi, transparent=True)
        if self.preview is not None:
            self.preview.mkdir(parents=True, exist_ok=True)
            fig.savefig(self.preview / f"{name}.png", dpi=PREVIEW_DPI,
                        facecolor="white")
        plt.close(fig)
        w, h = png_size(path)
        self.panels[name] = {"file": path.name,
                             "w": round(w / self.dpi * self.scale, 1),
                             "h": round(h / self.dpi * self.scale, 1)}
        logger.info("%-26s %4d x %4d px", name, w, h)

    def write_manifest(self, extra: dict | None = None) -> Path:
        """``panels.json`` next to the images: the sizes plus whatever the
        drawing needs to know about the data (legend entries, chosen K, where a
        heat map's rows sit)."""
        import json

        manifest = {"scale": self.scale, "panels": self.panels}
        manifest.update(extra or {})
        self.out.mkdir(parents=True, exist_ok=True)
        path = self.out / "panels.json"
        path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        logger.info("wrote %d panels to %s", len(self.panels), self.out)
        return path
