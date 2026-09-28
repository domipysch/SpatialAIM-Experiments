# Assemble the Results 2.2 figure as an Excalidraw drawing in the Obsidian vault.
#
# The data panels are the PNGs written by 02b_panels.py, embedded by wikilink --
# re-running that script updates them in place. This script draws what makes the
# panels readable as a pipeline: the arrows between them, the label each arrow
# carries, the panel letters and three legends. Anything else belongs in the
# caption.
#
# One mapper per invocation (``--mapper``, required), assembled from
# figure_2/<mapper>/panels/ into figure_2/<mapper>/figure_2_<mapper>.excalidraw.md.
#
# A band of keys, then three rows of four even columns, one arrow per column:
#
#   Induced type split  * ...      Criteria  - ...  - ...
#   * ...  * ...                   - ...  - ...
#
#   a               b                c                 d
#   | mapper        | mapper         | SpatialAIM w. mapper   | SpatialAIM w. mapper
#   v               v                v                 v
#   e               f                g                 h
#
#   i ...           the composition panel     interpreted cell types
#
# The two keys to the panels open the figure, side by side and two lines each --
# stacked they were seven entries deep and needed a column of the last row to
# themselves. The first row is what goes in -- the reference under each
# annotation, and the criterion scores of each SpatialAIM run; the second is the
# slice it comes out as. What each run was handed is in the panel headings, so
# the arrows only need to say which mapper made the step. The last row is the
# composition of the states in the fine annotation, named on its left and keyed
# on its right.
#
# Element construction, text measurement and the file format itself live in
# spatialaim_experiments.excalidraw; what is here is this figure's own layout.

import argparse
import json
import logging

from spatialaim_experiments import paths
from spatialaim_experiments.excalidraw import WIDTH, Drawing, text_size, wrap
from spatialaim_experiments.style import GREY, INK, MAPPERS, SCALE

logger = logging.getLogger(__name__)

FIGURE = "2"

PAD = 14
# The paper places the figure at its 5.17 in text width, so the drawing is scaled
# by 0.71 and every size is authored 1/0.71 larger than it should print: the
# letters come out at 8.6 pt on the page, everything else at 6.7 pt.
FS_LETTER = 50
FS_LEG = 40
LEG_PITCH = 48        # a row of legend is 48 tall
LEG_INDENT = 54       # swatch to label
LEG_GAP = 110         # one entry of a legend line to the next
LEG_COL_GAP = 80      # the two keys of the top band to each other
LEG_BAND_GAP = 46     # the band of keys at the top to the first row of letters
SWATCH = 22
LETTER_DROP = FS_LETTER * 1.2 + 6   # a letter sits this far above its panel

FS_HEAD = 40          # the heading next to a panel letter
HEAD_DX = 44          # letter to heading
COL_W = 470           # four even columns: two aligner baselines, two SpatialAIM runs
ARROW_LEN = 88        # the arrows between the first two rows
ROW3_GAP = 56         # second row's panels to the last row. Small, because
                      # nothing stands over that row any more: its letter and
                      # heading sit beside the panel, not above it.
INSET = 16            # an arrow starts this far below the panel it leaves
ARROW_CLEAR = 44      # and stops this far above the heading of the row below
SCORE_DROP = 6        # the score panels hang a little below the row's line, as
                      # placed by hand: their axes then sit level with the UMAPs

LINE = FS_HEAD * 1.2  # one line of heading

# What each panel says -- the wording is the figure's, so it lives here and not in
# the caption. A heading is centred on its letter: a one-line heading sits next
# to it, a two-line one straddles it half a line up. The letters of a row are all
# on the letter line either way.
HEADINGS = {
    "umap_celltype": "Coarse annotation",
    "umap_split": "Induced split",
    "scores_split": "SpatialAIM on split\ncell types",
    "scores_overcluster": "SpatialAIM without\nstart annotation",
    "map_celltype": "Mapped cell types",
    "map_split": "Glut. A & B mapped at random",
    "map_split_aim": "Glut. A & B re-merged",
    "map_overcluster_aim": "Finer cell types detected",
    "dots_states_fine": "Found types correspond to finer cell types",
}

# The two keys to the panels above. They name a band across the top of the
# figure now, not a column of the last row. The third belongs to (i) and stays
# beside it, set the same way.
HEAD_SPLIT = "Induced type split"
HEAD_CRITERIA = "Criteria"
HEAD_STATES = "interpreted cell types"

drawing: Drawing | None = None  # set by main(), once --mapper says which figure this is


def text(body: str, x, y, size=FS_LEG, colour=INK, right=False):
    return drawing.text(body, x, y, size=size, colour=colour,
                        align="right" if right else "left")


def swatch(x: float, y: float, entry, kind: str) -> None:
    """The mark in front of a legend label: a dot for a cell type, a rule in the
    curve's own weight for a criterion."""
    if kind == "dot":
        drawing.dot(x + SWATCH / 2, y + FS_LEG * 0.6, SWATCH, entry[1])
    else:
        drawing.line(x, y + FS_LEG * 0.6, 44, entry[1],
                     stroke=1.6 if entry[2] > 1 else 1)


def entry_width(entry) -> float:
    return LEG_INDENT + text_size(entry[0], FS_LEG)[0]


def legend_natural(heading: str, entries) -> float:
    """What a key would span if it were set on one line -- what its column's
    share of the width is worked out from."""
    return (text_size(heading, FS_LEG)[0] + LEG_GAP
            + sum(entry_width(entry) for entry in entries)
            + LEG_GAP * (len(entries) - 1))


def legend_block(x: float, y: float, width: float, heading: str, entries,
                 kind: str) -> float:
    """One key as a column: the heading, then its entries flowed across
    ``width`` and wrapped onto as many lines as they need. Returns the y below
    it.

    Stacked one entry per row, these two keys were seven rows deep and needed a
    column of the last row to themselves. Flowed, each is two lines, which is
    what lets them sit side by side above the figure. A wrapped line starts at
    the column's own left edge rather than under the first entry: the headings
    are wide, and indenting to them costs a line."""
    text(heading, x, y, colour=GREY)
    cursor = x + text_size(heading, FS_LEG)[0] + LEG_GAP
    line = y
    for entry in entries:
        w = entry_width(entry)
        if cursor > x and cursor + w > x + width:
            line, cursor = line + LEG_PITCH, x
        swatch(cursor, line, entry, kind)
        text(entry[0], cursor + LEG_INDENT, line)
        cursor += w + LEG_GAP
    return line + LEG_PITCH


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
    if manifest["k_split"] >= 4:
        HEADINGS["map_split_aim"] = "Split halves kept"

    def size(name: str) -> tuple[float, float]:
        return panels[name]["w"], panels[name]["h"]

    # even columns: the gaps are what is left over, and the two that carry an
    # SpatialAIM arrow have to hold its label as well
    gap = (WIDTH - 2 * PAD - 4 * COL_W) / 3
    left = [PAD + i * (COL_W + gap) for i in range(4)]
    middle = [x + COL_W / 2 for x in left]

    def place(name: str, column: int, top: float, letter: str,
              drop: float = 0, room: float = COL_W) -> dict:
        heading = wrap(HEADINGS[name], room - HEAD_DX, FS_HEAD)
        rise = heading.count("\n") * LINE / 2
        w, _ = size(name)
        x = middle[column] - w / 2
        element = drawing.image(name, x, top + drop, manifest)
        text(letter, x, top - LETTER_DROP, size=FS_LETTER)
        text(heading, x + HEAD_DX, top - LETTER_DROP + 9 - rise, size=FS_HEAD)
        return element

    # A heading reaches this far above its letter line -- half a line per extra
    # line of it, since the heading is centred on the letter. It is what the space
    # over a row has to hold, and where an arrow into that row has to stop.
    def heading_lift(name: str) -> float:
        return wrap(HEADINGS[name], COL_W - HEAD_DX, FS_HEAD).count("\n") * LINE / 2

    # ---------------- the two keys, as a band across the top ----------------
    # They used to stand in the last row, seven stacked entries deep. Side by
    # side, each flowed over two lines, they cost a quarter of that height and
    # can sit above the figure, where they are read before the panels rather
    # than after them. Each column gets the share of the width it would need on
    # one line, so the two wrap equally hard.
    keys = [(HEAD_SPLIT, legends["annotation"], "dot"),
            (HEAD_CRITERIA, legends["criteria"], "line")]
    want = [legend_natural(heading, entries) for heading, entries, _ in keys]
    room = WIDTH - 2 * PAD - LEG_COL_GAP * (len(keys) - 1)
    band, cursor = PAD, PAD
    for (heading, entries, kind), share in zip(keys, want):
        width = room * share / sum(want)
        band = max(band, legend_block(cursor, PAD, width, heading, entries, kind))
        cursor += width + LEG_COL_GAP

    # The score panels are a fraction of a unit taller than the UMAPs, so the row
    # ends where the deepest of its panels does and every arrow starts there. Each
    # row is preceded by the same free space plus its band of letters, so the two
    # gaps down the figure are equal.
    row1_names = ("umap_celltype", "umap_split", "scores_split",
                  "scores_overcluster")
    row1_top = (band + LEG_BAND_GAP + LETTER_DROP - 9
                + max(heading_lift(name) for name in row1_names))
    row1_bottom = row1_top + max(size("umap_celltype")[1], size("umap_split")[1],
                                 size("scores_split")[1] + SCORE_DROP,
                                 size("scores_overcluster")[1] + SCORE_DROP)

    row2_names = ("map_celltype", "map_split", "map_split_aim",
                  "map_overcluster_aim")
    row2_top = (row1_bottom + INSET + ARROW_LEN + ARROW_CLEAR + LETTER_DROP - 9
                + max(heading_lift(name) for name in row2_names))
    row2_bottom = row2_top + size("map_celltype")[1]
    row3_top = row2_bottom + ROW3_GAP

    # The letters run row by row, a-b-c-d over e-f-g-h: what goes in, then what
    # comes out of it.
    # ---------------- first row: the reference, and what SpatialAIM scores ----------
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
    name = manifest["mapper_name"]
    for column in range(4):
        # the label always sits to the left of its arrow, and the two of them
        # together are centred in the column -- so the arrow itself is off centre
        # by half the label
        # "SpatialAIM with Tangram" is wider than its column but the columns are
        # not touching, so it is allowed to reach into the gap and stay on one
        # line, as set by hand in the drawing. A long mapper name ("SpatialAIM
        # with nearest centroid") is wider than even that and still breaks over
        # two lines rather than running into the neighbouring panel.
        room = (COL_W if column < 2 else COL_W + gap) - 60
        label = wrap(name if column < 2 else f"SpatialAIM with {name}", room, FS_LEG)
        x = middle[column] + (text_size(label, FS_LEG)[0] + 12) / 2
        end = (row2_top - LETTER_DROP + 9 - heading_lift(row2_names[column])
               - ARROW_CLEAR)
        drawing.arrow([(x, start), (x, end)])
        text(label, x - 12, start + 10, right=True)

    # ---------------- last row: the composition panel and its key -----------
    # Three blocks across the row, as placed by hand: the panel's letter and
    # heading in the first column, the panel itself in the middle, its key where
    # the last column of the rows above begins -- so the row is read against the
    # same grid even though it is not built from it.
    states = [(f"{s['n']}  {s['name']}", s["colour"]) for s in legends["states"]]
    head = HEADINGS["dots_states_fine"]
    key_x = left[3]
    head_right = PAD + HEAD_DX + text_size(head, FS_HEAD)[0]
    dots = drawing.image("dots_states_fine",
                 head_right + (key_x - head_right - size("dots_states_fine")[0]) / 2,
                 row3_top, manifest)

    # The state key is not a stack of rows but a key to the panel: each entry
    # sits at the height of the row it names, and its heading -- one line, which
    # this column has room for -- stands over the first of them.
    row3_h = size("dots_states_fine")[1]
    first = row3_top + legends["states"][0]["y"] * row3_h - FS_LEG * 0.6
    heading_y = first - text_size(HEAD_STATES, FS_LEG)[1] - LINE / 2
    text(HEAD_STATES, key_x, heading_y, colour=GREY)
    for entry, state in zip(states, legends["states"]):
        y = row3_top + state["y"] * row3_h - FS_LEG * 0.6
        drawing.dot(key_x + SWATCH / 2, y + FS_LEG * 0.6, SWATCH, entry[1])
        text(entry[0], key_x + LEG_INDENT, y)

    # (i) is named beside its panel rather than over it, on the line its key
    # opens on.
    text("i", PAD, heading_y, size=FS_LETTER)
    text(head, PAD + HEAD_DX, heading_y + 9, size=FS_HEAD)
    logger.info("canvas %.0f x %.0f units", WIDTH, dots["y"] + dots["height"] + PAD)


def main() -> None:
    global drawing
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mapper", required=True, choices=MAPPERS,
                    help="which mapper's panels to assemble")
    mapper = ap.parse_args().mapper
    panel_dir = paths.panel_dir(FIGURE, mapper)
    drawing = Drawing("f22", panel_dir, paths.VAULT)

    build(json.loads((panel_dir / "panels.json").read_text(encoding="utf-8")))
    drawing.check()
    drawing.write(paths.drawing(FIGURE, mapper))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
