# Assemble the Results 2.3 reconstruction figure as an Excalidraw drawing in the
# Obsidian vault.
#
# Same construction as figure_2_2: the data panels are the PNGs written by
# panels.py and embedded by wikilink, so re-running that script updates them in
# place. This script draws only what makes them readable -- panel letters,
# headings and the legends. Anything else belongs in the caption.
#
# Two rows:
#
#   a  AIM vs annotation, per spot | per gene | cell types used
#   b  the five mappers at the K each proposes
#
# The first row is one argument in three panels of identical geometry:
# reconstruction barely moves either way, and AIM gets there with substantially
# fewer cell types. The second asks which mapper does it best.
#
# Element construction, text measurement and the file format itself live in
# aim_experiments.excalidraw; what is here is this figure's own layout.

import json
import logging

from aim_experiments import paths
from aim_experiments.excalidraw import WIDTH, Drawing, advance, text_size

logger = logging.getLogger(__name__)

FIGURE = "2_3"
PANELS = paths.panel_dir(FIGURE)
OUT = paths.drawing(FIGURE)

PAD = 14
# The paper places the figure at its 5.17 in text width, so the drawing is
# scaled by 0.71 and every size is authored 1/0.71 larger than it should print.
FS_LETTER = 50
FS_HEAD = 40
FS_LEG = 40
FS_ROW = 30           # the per-row linkage labels, smaller than the mapper name
ROW_LABEL_W = sum(advance(c) for c in "ward") * FS_ROW   # "ward", the longest of them, right-aligned
LEG_PITCH = 48
# The mapper swatch, in both rows: a rounded bar tall enough to span most of the
# two rows one mapper occupies in row (a). Those rows are 76 units together, so
# 48 covers them without reaching into the neighbouring mapper's pair, and it is
# also exactly one line of FS_LEG, so it sits flush with its label in the key.
MARK_W, MARK_H = 16, 48
LEG_INDENT = MARK_W + 18
HEAD_DX = 44
LETTER_DROP = FS_LETTER * 1.2 + 6
LINE = FS_HEAD * 1.2

ROW_TOP = 50
ROW_GAP = 52          # panel bottom to the next row's letter line
COL_GAP = 60

# One heading per row. The three heat maps say what they measure on their own
# colour bars, so the row only has to say what the comparison is.
HEADINGS = {
    "delta_spot": "Reconstruction and granularity, comparing AIM and given annotations",
    "mappers": "Comparing different mappers within AIM",
}

drawing = Drawing("f23", PANELS, paths.VAULT)


def text(body: str, x, y, size=FS_LEG, **kw):
    return drawing.text(body, x, y, size=size, **kw)


def bar(cx, cy, w, h, colour):
    """The mapper swatch of both rows: a rounded bar, centred on (cx, cy).
    Taller than the dot it replaced, so in row (a) it reads as spanning the pair
    of rows it names rather than pointing at the gap between them."""
    return drawing.bar(cx, cy, w, h, colour)


def legend_row(x: float, y: float, entries, gap: float = 34) -> float:
    """The same key laid out left to right. The two heat maps together leave no
    column free beside them, so a stacked key has nowhere to go.

    Its swatches are the same rounded bars row (a) marks its mappers with, at the
    same size, so the two are read as one another without a second look. A bar of
    ``MARK_H`` is exactly one line high, so it sits flush with its label."""
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
    legends = manifest["legends"]

    def size(name: str) -> tuple[float, float]:
        return panels[name]["w"], panels[name]["h"]

    def label(name: str, x: float, top: float, letter: str) -> None:
        """The row's letter and its single heading, on one line at x."""
        text(letter, x, top - LETTER_DROP, size=FS_LETTER)
        text(HEADINGS[name], x + HEAD_DX, top - LETTER_DROP + 9, size=FS_HEAD)

    # Both rows occupy the same span: they start at PAD and end where the wider
    # of the two does. Row (b) is one image, so its width sets the span, and
    # row (a)'s mapper-name column takes up whatever the three heat maps leave.
    w_map, h_map = size("mappers")
    w_heat, h_heat = size("delta_spot")
    gap = 30
    span = w_map
    names_w = span - 3 * w_heat - 2 * gap    # the column of mapper names
    heat_x = PAD + names_w

    # ---------------- row 1: the corpus --------------------------------------
    # Three heat maps of identical geometry. Their columns are labelled inside
    # the images and their colour bars name their own units, so the row needs
    # one heading and nothing else. The mapper names are drawn once, here, to
    # the left of the row: inside the first image they would have made its cells
    # narrower than the other two's.
    row1 = ROW_TOP + LETTER_DROP
    label("delta_spot", PAD, row1, "a")
    for i, name in enumerate(("delta_spot", "delta_gene", "delta_k")):
        drawing.image(name, heat_x + i * (w_heat + gap), row1, manifest)

    # Each mapper occupies two rows -- one per linkage -- so the column carries
    # both: the linkage of every row against the heat map, and, centred on the
    # pair, a dot in the mapper's colour. The dot rather than the name, because
    # row (b)'s key already spells the five out and repeating them here cost more
    # width than the heat maps could spare.
    for tree, y in zip(manifest["trees"] * 5, manifest["heat_row_lines"]):
        w, h = text_size(tree, FS_ROW)
        text(tree, heat_x - 12 - w, row1 + y * h_heat - h / 2, size=FS_ROW)

    mark_x = heat_x - 12 - ROW_LABEL_W - 16 - MARK_W / 2
    for entry, y in zip(legends["mappers"], manifest["heat_rows"]):
        bar(mark_x, row1 + y * h_heat, MARK_W, MARK_H, entry["colour"])

    # ---------------- row 2: the mappers -------------------------------------
    row2 = row1 + h_heat + ROW_GAP + LETTER_DROP
    label("mappers", PAD, row2, "b")
    drawing.image("mappers", PAD, row2, manifest)

    # the only key the figure needs: which colour is which mapper
    y = row2 + h_map + 24
    legend_row(PAD, y, [(e["label"], e["colour"]) for e in legends["mappers"]])
    logger.info("canvas %.0f x %.0f units", WIDTH, y + LEG_PITCH + PAD)


def main() -> None:
    manifest = json.loads((PANELS / "panels.json").read_text(encoding="utf-8"))
    build(manifest)
    drawing.check()
    drawing.write(OUT)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
