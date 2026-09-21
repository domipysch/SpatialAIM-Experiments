# Assemble the Results 2.2 figure as an Excalidraw drawing in the Obsidian vault.
#
# The data panels are the PNGs written by panels.py, embedded by wikilink --
# re-running that script updates them in place. This script draws what makes the
# panels readable as a pipeline: the arrows between them, the label each arrow
# carries, the panel letters and three legends. Anything else belongs in the
# caption.
#
# One mapper per invocation (``--mapper``, default Tangram, the figure the paper
# prints); every other mapper is assembled from figure_2_2_<mapper>/ into
# figure_2_2_<mapper>.excalidraw.md.
#
# Three rows, four even columns, one arrow per column:
#
#   a               b                c                 d
#   | mapper        | mapper         | AIM w. mapper   | AIM w. mapper
#   v               v                v                 v
#   e               f                g                 h
#
#   cell types      criteria         i                 cell states
#
# The first row is what goes in -- the reference under each annotation, and the
# criterion scores of each AIM run; the second is the slice it comes out as. What
# each run was handed is in the panel headings, so the arrows only need to say
# which mapper made the step. The last row is the legends and, between them, the
# composition of the seven states in the fine annotation.
#
# Element construction, text measurement and the file format itself live in
# aim_experiments.excalidraw; what is here is this figure's own layout.

import argparse
import json
import logging

from aim_experiments import paths
from aim_experiments.excalidraw import WIDTH, Drawing, text_size, wrap
from aim_experiments.style import GREY, INK, MAPPERS, SCALE

logger = logging.getLogger(__name__)

FIGURE = "2_2"
# set by main() from --mapper; Tangram keeps the unsuffixed names, every other
# mapper gets its own copy of the same figure next to it
PANELS = paths.panel_dir(FIGURE)
OUT = paths.drawing(FIGURE)

PAD = 14
# The paper places the figure at its 5.17 in text width, so the drawing is scaled
# by 0.71 and every size is authored 1/0.71 larger than it should print: the
# letters come out at 8.6 pt on the page, everything else at 6.7 pt.
FS_LETTER = 50
FS_LEG = 40
LEG_PITCH = 48        # a row of legend is 48 tall
LEG_INDENT = 54       # swatch to label
SWATCH = 22
LETTER_DROP = FS_LETTER * 1.2 + 6   # a letter sits this far above its panel

FS_HEAD = 40          # the heading next to a panel letter
HEAD_DX = 44          # letter to heading
COL_W = 470           # four even columns: two aligner baselines, two AIM runs
ROW1_TOP = 50
ARROW_LEN = 88        # the arrows between the first two rows
ROW3_GAP = 130        # second row's panels to the last row
INSET = 16            # an arrow starts this far below the panel it leaves
ARROW_CLEAR = 22      # and stops this far above the heading of the row below
SCORE_DROP = 6        # the score panels hang a little below the row's line, as
                      # placed by hand: their axes then sit level with the UMAPs

LINE = FS_HEAD * 1.2  # one line of heading

# What each panel says -- the wording is the figure's, so it lives here and not in
# the caption. The second value is how far the heading is lifted off the letter
# line; None, as they all are while every heading is a single line, puts its last
# line on that line, next to the letter.
HEADINGS = {
    "umap_celltype": ("Coarse cell types", None),
    "umap_split": ("Invisible split", None),
    "scores_split": ("AIM on split" + chr(10) + "cell types", None),
    "scores_overcluster": ("AIM without" + chr(10) + "start annotation", None),
    "map_celltype": ("Mapped cell types", None),
    "map_split": ("Split halves mapped at random", None),
    "map_split_aim": ("Split halves re-merged", None),
    "map_overcluster_aim": ("Finer cell types detected", None),
    "dots_states_fine": ("Found types correspond to finer cell types", None),
}

# set by main(), once --mapper says which figure this is
drawing = Drawing("f22", PANELS, paths.VAULT)


def dot(cx, cy, d, colour):
    return drawing.dot(cx, cy, d, colour)


def text(body: str, x, y, size=FS_LEG, colour=INK, right=False):
    return drawing.text(body, x, y, size=size, colour=colour,
                        align="right" if right else "left")


def image(name: str, x, y, manifest: dict):
    return drawing.image(name, x, y, manifest)


def arrow(points, **kw):
    return drawing.arrow(points, **kw)



def legend_width(entries) -> float:
    return LEG_INDENT + max(text_size(entry[0], FS_LEG)[0] for entry in entries)


def legend(x: float, y: float, entries, kind="dot") -> float:
    """Swatch-and-label rows; returns the y below the last row."""
    for entry in entries:
        label, colour = entry[0], entry[1]
        if kind == "dot":
            dot(x + SWATCH / 2, y + FS_LEG * 0.6, SWATCH, colour)
        else:
            drawing.line(x, y + FS_LEG * 0.6, 44, colour,
                         stroke=1.6 if entry[2] > 1 else 1)
        text(label, x + LEG_INDENT, y)
        y += LEG_PITCH
    return y


# --------------------------------------------------------------------------
# the figure
# --------------------------------------------------------------------------

def build(manifest: dict) -> None:
    panels = manifest["panels"]
    legends = manifest["legends"]

    # (g) is the level the split run proposes. That is the re-merge for every
    # mapper but DOT, which proposes the annotation it was handed (K = 4, the
    # four annotated types) -- there the panel shows the split kept, and saying
    # "re-merged" over it would contradict the picture underneath.
    if manifest.get("k_split", 3) >= 4:
        HEADINGS["map_split_aim"] = ("Split halves kept", None)

    def size(name: str) -> tuple[float, float]:
        return panels[name]["w"], panels[name]["h"]

    # even columns: the gaps are what is left over, and the two that carry an
    # AIM arrow have to hold its label as well
    gap = (WIDTH - 2 * PAD - 4 * COL_W) / 3
    left = [PAD + i * (COL_W + gap) for i in range(4)]
    middle = [x + COL_W / 2 for x in left]

    def place(name: str, column: int, top: float, letter: str,
              drop: float = 0, room: float = COL_W) -> dict:
        heading, rise = HEADINGS[name]
        heading = wrap(heading, room - HEAD_DX, FS_HEAD)
        if rise is None:
            rise = heading.count(chr(10)) * LINE
        w, _ = size(name)
        x = middle[column] - w / 2
        element = image(name, x, top + drop, manifest)
        text(letter, x, top - LETTER_DROP, size=FS_LETTER)
        text(heading, x + HEAD_DX, top - LETTER_DROP + 9 - rise, size=FS_HEAD)
        return element

    # The score panels are a fraction of a unit taller than the UMAPs, so the row
    # ends where the deepest of its panels does and every arrow starts there. Each
    # row is preceded by the same free space plus its band of letters, so the two
    # gaps down the figure are equal.
    row1_top = ROW1_TOP
    row1_bottom = row1_top + max(size("umap_celltype")[1], size("umap_split")[1],
                                 size("scores_split")[1] + SCORE_DROP,
                                 size("scores_overcluster")[1] + SCORE_DROP)
    # An arrow stops above the *heading* of the panel it points at, not above its
    # letter, so the gap between the rows has to hold the arrow, its two
    # clearances and however many lines those headings run to.
    def heading_lift(name: str) -> float:
        return wrap(HEADINGS[name][0], COL_W - HEAD_DX,
                    FS_HEAD).count(chr(10)) * LINE

    row2_names = ("map_celltype", "map_split", "map_split_aim",
                  "map_overcluster_aim")
    row2_top = (row1_bottom + INSET + ARROW_LEN + ARROW_CLEAR + LETTER_DROP - 9
                + max(heading_lift(name) for name in row2_names))
    row2_bottom = row2_top + size("map_celltype")[1]
    row3_top = row2_bottom + ROW3_GAP

    # The letters run row by row, a-b-c-d over e-f-g-h: what goes in, then what
    # comes out of it.
    # ---------------- first row: the reference, and what AIM scores ----------
    place("umap_celltype", 0, row1_top, "a")
    place("umap_split", 1, row1_top, "b")
    place("scores_split", 2, row1_top, "c", drop=SCORE_DROP)
    place("scores_overcluster", 3, row1_top, "d", drop=SCORE_DROP)

    # ---------------- second row: the slice, mapped four ways ---------------
    place("map_celltype", 0, row2_top, "e")
    place("map_split", 1, row2_top, "f")
    place("map_split_aim", 2, row2_top, "g")
    place("map_overcluster_aim", 3, row2_top, "h")

    # ---------------- the arrows ----------------
    # One per column, from what goes in to the slice it produces: all four run
    # down the middle of their column, over the same span, with the mapper's name
    # hung off the same side at the same offset.
    start = row1_bottom + INSET
    for column in range(4):
        # the label always sits to the left of its arrow, and the two of them
        # together are centred in the column -- so the arrow itself is off centre
        # by half the label
        name = manifest.get("mapper_name", "Tangram")
        # a long mapper name ("AIM with nearest centroid") is wider than the
        # column, so it breaks over two lines rather than into its neighbour
        label = wrap(name if column < 2 else f"AIM with {name}", COL_W - 60, FS_LEG)
        x = middle[column] + (text_size(label, FS_LEG)[0] + 12) / 2
        end = (row2_top - LETTER_DROP + 9 - heading_lift(row2_names[column])
               - ARROW_CLEAR)
        arrow([(x, start), (x, end)])
        text(label, x - 12, start + 10, right=True)

    # ---------------- last row: the legends and the composition panel -------
    annotation = legends["annotation"]
    criteria = legends["criteria"]
    states = [(f"{s['n']}  {s['name']}", s["colour"]) for s in legends["states"]]
    # at this text size four blocks no longer fit across the row, so the two keys
    # to the panels above are stacked into one
    blocks = [max(legend_width(annotation), legend_width(criteria)),
              size("dots_states_fine")[0], legend_width(states)]
    gap = (WIDTH - 2 * PAD - sum(blocks)) / (len(blocks) - 1)
    x, cursor = [], PAD
    for block in blocks:
        x.append(cursor)
        cursor += block + gap

    row3_h = size("dots_states_fine")[1]

    heading_y = row3_top - LETTER_DROP + 9
    y = heading_y
    for entries, heading, kind in [(annotation, "cell types", "dot"),
                                   (criteria, "criteria", "line")]:
        text(heading, x[0], y, colour=GREY)
        y = legend(x[0], y + LINE, entries, kind=kind) + LINE / 2

    dots = image("dots_states_fine", x[1], row3_top, manifest)
    text("i", x[1], row3_top - LETTER_DROP, size=FS_LETTER)
    text(HEADINGS["dots_states_fine"][0], x[1] + HEAD_DX,
         row3_top - LETTER_DROP + 9, size=FS_HEAD)

    # The state legend is not a stack of rows but a key to the panel: each entry
    # sits at the height of the row it names. Its heading is wider than the block,
    # so it does not share the heading line of the row -- it sits directly over
    # its own first entry, where it was placed by hand.
    first = row3_top + legends["states"][0]["y"] * row3_h - FS_LEG * 0.6
    # at this text size the heading no longer fits the block on one line, and it
    # would run into the panel's own bracket labels if it started further left
    heading = wrap("interpreted cell types", blocks[2], FS_LEG)
    text(heading, x[2],
         first - text_size(heading, FS_LEG)[1] - LINE / 2, colour=GREY)
    for entry, state in zip(states, legends["states"]):
        y = row3_top + state["y"] * row3_h - FS_LEG * 0.6
        dot(x[2] + SWATCH / 2, y + FS_LEG * 0.6, SWATCH, entry[1])
        text(entry[0], x[2] + LEG_INDENT, y)
    logger.info("canvas %.0f x %.0f units", WIDTH, dots["y"] + dots["height"] + PAD)


def main() -> None:
    global PANELS, OUT, drawing
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mapper", default="tangram", choices=MAPPERS,
                    help="which mapper's panels to assemble (default tangram)")
    mapper = ap.parse_args().mapper
    copy = None if mapper == "tangram" else mapper
    PANELS = paths.panel_dir(FIGURE, copy)
    OUT = paths.drawing(FIGURE, copy)
    drawing = Drawing("f22", PANELS, paths.VAULT)

    manifest = json.loads((PANELS / "panels.json").read_text(encoding="utf-8"))
    build(manifest)
    drawing.check()
    drawing.write(OUT)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
