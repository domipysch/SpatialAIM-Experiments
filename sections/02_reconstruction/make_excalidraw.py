# Assemble the Results 2.2 reconstruction figure as an Excalidraw drawing in the
# Obsidian vault.
#
# Same construction as figure 2: the data panels are the PNGs written by
# panels.py and embedded by wikilink, so re-running that script updates them in
# place. This script draws only what makes them readable -- panel letters,
# headings and the legends. Anything else belongs in the caption.
#
# The key to the five mapper colours, and below it three rows:
#
#   a  SpatialAIM vs annotation, per spot | per gene | cell types used
#   b  the five mappers at the K each proposes, per cell   c per gene
#      d runtime   e memory
#   f  the start clusters   g ward cut   h average cut   (UMAPs)
#
# The first row is one argument in three panels of identical geometry:
# reconstruction barely moves either way, and SpatialAIM gets there with substantially
# fewer cell types. The second asks which mapper does it best. The third shows
# on one reference what the two linkages make of the same start clusters.
#
# The UMAPs of the third row are screenshots taken by hand from the GUI, not
# written by panels.py: they live in figure_2/umaps/ and are placed at the row's
# width, at their own aspect ratio.
#
# Element construction, text measurement and the file format itself live in
# spatialaim_experiments.excalidraw; what is here is this figure's own layout.

import json
import logging

from spatialaim_experiments import paths
from spatialaim_experiments.excalidraw import WIDTH, Drawing, advance, text_size
from spatialaim_experiments.style import (HEAT_AXES, MAPPER_COLOURS, MAPPER_LABELS, MAPPERS,
                                          TREE_LABELS, TREES, png_size)

logger = logging.getLogger(__name__)

FIGURE = "2"
PANELS = paths.panel_dir(FIGURE)
OUT = paths.drawing(FIGURE)
UMAPS = PANELS.parent / "umaps"   # the hand-taken screenshots of row (f-h)

PAD = 14
# The paper places the figure at its 5.17 in text width, so the drawing is
# scaled by 0.71 and every size is authored 1/0.71 larger than it should print.
FS_LETTER = 50
FS_HEAD = 40
FS_LEG = 40
FS_ROW = 30           # the per-row linkage labels, smaller than the mapper name
# The width of the column those labels are right-aligned in: the widest of them,
# "ward". text_size over-measures Helvetica slightly, so the box always holds the
# label and never wraps it.
ROW_LABEL_W = sum(advance(c) for c in "ward") * FS_ROW
# The mapper swatch: a rounded vertical bar. In row (a) it is drawn at the height
# of the two heat-map rows the mapper occupies, so the
# five bars tile the column edge to edge with no gaps between them. In the key it
# is one line of FS_LEG instead, which sits flush with its label.
MARK_W, MARK_H = 16, 48
LEG_INDENT = MARK_W + 18
HEAD_DX = 44
LETTER_DROP = FS_LETTER * 1.2 + 6
LINE = FS_HEAD * 1.2

ROW_TOP = 50
ROW_GAP = 28          # row (a)'s panel bottom to row (b-e)'s letter line
ROW3_GAP = 20         # the same from row (b-e) to row (f-h)
COL_GAP = 60
LEG_GAP = 30          # key baseline to the letter line of row (a) below it

# One heading per row. The three heat maps say what they measure on their own
# colour bars, so the row only has to say what the comparison is.
HEADINGS = {
    "delta_spot": "Reconstruction and granularity, comparing SpatialAIM and given annotations",
    "mappers_cell": "Comparing mappers, per cell",
    "mappers_gene": "Mappers, per gene",
    "runtime": "Runtime, all K",
    "memory": "Peak memory",
    "start": "Leiden start clusters",
    "ward": "Ward linkage, K = 2",
    "average": "Average linkage, K = 2",
}
UMAP_ROW = ("start", "ward", "average")

drawing = Drawing("f22", PANELS, paths.VAULT)


def text(body: str, x, y, size=FS_LEG, **kw):
    return drawing.text(body, x, y, size=size, **kw)


def bar(cx, cy, w, h, colour):
    """The mapper swatch: a rounded bar, centred on (cx, cy). In row (a) its
    height is the pair of heat-map rows it names, so it brackets them rather than
    pointing at the gap between them."""
    return drawing.bar(cx, cy, w, h, colour)


def legend_row(x: float, y: float, entries, gap: float = 34) -> float:
    """The key, laid out left to right above the figure: it names the colours
    both rows use, so it belongs ahead of them rather than under the last panel.
    The heat maps leave no column free beside them, so a stacked key has nowhere
    to go anyway.

    Its swatches are the same rounded bars row (a) marks its mappers with, so the
    two are read as one another without a second look -- shorter here, because a
    bar of ``MARK_H`` is exactly one line high and sits flush with its label."""
    cursor = x
    for label, colour in entries:
        bar(cursor + MARK_W / 2, y + FS_LEG * 0.6, MARK_W, MARK_H, colour)
        text(label, cursor + LEG_INDENT, y)
        cursor += LEG_INDENT + text_size(label, FS_LEG)[0] + gap
    return cursor - gap


# --------------------------------------------------------------------------
# the figure
# --------------------------------------------------------------------------

def build(manifest: dict) -> None:
    panels = manifest["panels"]

    def size(name: str) -> tuple[float, float]:
        return panels[name]["w"], panels[name]["h"]

    def label(name: str, x: float, top: float, letter: str) -> None:
        """The row's letter and its single heading, on one line at x."""
        text(letter, x, top - LETTER_DROP, size=FS_LETTER)
        text(HEADINGS[name], x + HEAD_DX, top - LETTER_DROP + 9, size=FS_HEAD)

    # All rows occupy the same span: they start at PAD and end where row (b-e) --
    # the mapper comparison and the two cost panels -- does; row (a)'s
    # mapper-name column takes up whatever the three heat maps leave.
    gap = 30
    w_cell, h_map = size("mappers_cell")
    w_map = w_cell + gap + size("mappers_gene")[0]   # b and c together
    w_cost = size("runtime")[0]
    w_heat, h_heat = size("delta_spot")
    span = w_map + 2 * w_cost + 2 * gap
    names_w = span - 3 * w_heat - 2 * gap    # the column of mapper names
    heat_x = PAD + names_w

    # ---------------- the key ------------------------------------------------
    # Which colour is which mapper, said once, at the top: both rows use the five
    # colours, so the reader needs them before panel (a), not after panel (b).
    legend_row(PAD, ROW_TOP, [(MAPPER_LABELS[m], MAPPER_COLOURS[m]) for m in MAPPERS])

    # ---------------- row 1: the corpus --------------------------------------
    # Three heat maps of identical geometry. Their columns are labelled inside
    # the images and their colour bars name their own units, so the row needs
    # one heading and nothing else. The mapper names are drawn once, here, to
    # the left of the row: inside the first image they would have made its cells
    # narrower than the other two's.
    row1 = ROW_TOP + FS_LEG * 1.2 + LEG_GAP + LETTER_DROP
    label("delta_spot", PAD, row1, "a")
    for i, name in enumerate(("delta_spot", "delta_gene", "delta_k")):
        drawing.image(name, heat_x + i * (w_heat + gap), row1, manifest)

    # Each mapper occupies two rows -- one per linkage -- so the column carries
    # both: the linkage of every row against the heat map, right-aligned on the
    # image so the two words end on one edge, and, spanning the pair exactly, a
    # bar in the mapper's colour. The bar rather than the name, because the key
    # already spells the five out and repeating them here cost more width than
    # the heat maps could spare.
    # The heat map's rows (mapper x tree, mapper by mapper) split HEAT_AXES evenly;
    # HEAT_AXES counts from the image's bottom, the canvas from its top (row1).
    _, bottom, _, height = HEAT_AXES
    row_h = height * h_heat / (len(MAPPERS) * len(TREES))
    rows_top = row1 + (1 - bottom - height) * h_heat
    for i in range(len(MAPPERS) * len(TREES)):
        tree = TREE_LABELS[TREES[i % len(TREES)]]
        h = text_size(tree, FS_ROW)[1]
        text(tree, heat_x - 12, rows_top + (i + 0.5) * row_h - h / 2, size=FS_ROW,
             align="right", box=ROW_LABEL_W)
    mark_x = heat_x - 12 - ROW_LABEL_W - 16 - MARK_W / 2
    pair_h = len(TREES) * row_h
    for i, mapper in enumerate(MAPPERS):
        bar(mark_x, rows_top + (i + 0.5) * pair_h, MARK_W, pair_h, MAPPER_COLOURS[mapper])

    # ---------------- row 2: the mappers -------------------------------------
    row2 = row1 + h_heat + ROW_GAP + LETTER_DROP
    for name, x, letter in (("mappers_cell", PAD, "b"),
                            ("mappers_gene", PAD + w_cell + gap, "c")):
        label(name, x, row2, letter)
        drawing.image(name, x, row2, manifest)
    for i, (name, letter) in enumerate((("runtime", "d"), ("memory", "e"))):
        x = PAD + w_map + gap + i * (w_cost + gap)
        label(name, x, row2, letter)
        drawing.image(name, x, row2, manifest)

    # ---------------- row 3: the UMAPs ---------------------------------------
    # Three screenshots across the same span as the rows above, with the same
    # gap, each under its own letter and heading.
    row3 = row2 + h_map + ROW3_GAP + LETTER_DROP
    add_umaps(manifest, (span - 2 * gap) / 3)
    h_umap = 0
    for i, (name, letter) in enumerate(zip(UMAP_ROW, "fgh")):
        x = PAD + i * (panels[name]["w"] + gap)
        label(name, x, row3, letter)
        drawing.image(name, x, row3, manifest, folder=UMAPS)
        h_umap = max(h_umap, panels[name]["h"])
    logger.info("canvas %.0f x %.0f units", WIDTH, row3 + h_umap + PAD)


def add_umaps(manifest: dict, width: float) -> None:
    """Enter the UMAP screenshots into the manifest as if panels.py had written
    them, scaled to ``width`` at their own aspect ratio."""
    for name in UMAP_ROW:
        w, h = png_size(UMAPS / f"{name}.png")
        manifest["panels"][name] = {"file": f"{name}.png",
                                    "w": width, "h": round(width * h / w, 1)}


def main() -> None:
    manifest = json.loads((PANELS / "panels.json").read_text(encoding="utf-8"))
    build(manifest)
    drawing.check()
    drawing.write(OUT)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
