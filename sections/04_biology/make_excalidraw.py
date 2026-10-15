# Assemble the Results 2.4 figure -- biological and topological plausibility --
# as an Excalidraw drawing in the Obsidian vault.
#
# Same construction as figures 2 and 3: the data panels are the PNGs
# written by panels.py and embedded by wikilink, so re-running that script
# updates them in place. This script draws only what makes them readable --
# letters, headings and the two words naming each sketch half. Anything else
# belongs in the caption.
#
# One mapper per invocation (``--mapper``; all mappers are equal), or the version
# aggregated over all five (``--all-mappers``):
# figure_4/<m>/panels/ -> figure_4/<m>/figure_4_<m>.excalidraw.md, m = mapper or all_mappers.
#
# An optional top row, then two rows of three, one row per criterion, and every
# cell its own panel:
#
#   a kidney, annotation | b kidney, SpatialAIM | c wins
#   d sketch | e neighbourhood enrichment | f local spatial purity
#   g sketch | h modularity over the slice | i modularity over the reference
#
# The top row is the worked example and the overview: one tissue mapped both
# ways (04_kidney_example/, whose panels.json also supplies each map's key and
# metrics), then on how many references SpatialAIM beats the annotation (one bar
# per metric). It is drawn in one version only: by default nearest centroid's
# (``--top-row`` / ``--no-top-row`` override). Without it the letters start at
# the sketch, a-f. Letters are handed out in reading order, so they never need
# renumbering by hand.
#
# Each row opens with a sketch of what its criterion means, drawn on a toy graph
# with the same labels twice -- observed above, permuted below. Then the two
# measures of that criterion, SpatialAIM against the author's annotation, the diagonal
# being parity. The single colour bar under the block applies to all four.
#
# The letter and its caption sit at the top left of the cell they belong to. The
# three letters of a row share one line, and the images share the line below the
# row's tallest caption, so both read as a row however long a caption is.
#
# Run with any Python (no scientific dependencies).

# Element construction, text measurement and the file format itself live in
# spatialaim_experiments.excalidraw; what is here is this figure's own layout.

import argparse
import json
import logging

from spatialaim_experiments import paths
from spatialaim_experiments.excalidraw import WIDTH, Drawing, text_size, wrap
from spatialaim_experiments.style import GREY, MAPPERS, SCALE

logger = logging.getLogger(__name__)

FIGURE = "4"

PAD = 14
# The paper places the figure at its 5.17 in text width, so the drawing is
# scaled by 0.71 and every size is authored 1/0.71 larger than it should print.
FS_LETTER = 50
FS_HEAD = 40
HEAD_DX = 44
LINE = FS_HEAD * 1.2
# The caption's baseline sits this far above the image it belongs to.
HEAD_GAP = 6

ROW_TOP = 10
ROW_GAP = 40          # panel bottom to the next row's letter line
BAR_GAP = 26

# One letter and one caption per cell. The two sketches keep the caption that
# used to head their whole row -- they are that criterion made visual -- and the
# four scatters keep the heading that used to stand over them.
#
# The wording is the user's own, read back out of the drawing after they edited
# it there; keep it verbatim rather than paraphrasing on the next regeneration.
# So are the line breaks, which balance the two lines where a greedy wrap would
# leave one word alone -- ``wrap`` keeps a break that is already in the text and
# is only the fallback for a caption that still does not fit its cell.
TOP_ROW_DEFAULT = {"nearest_centroid"}   # the versions that get the top row
# The top row's captions; {k} is the map's number of cell types.
HEAD_WINS = "SpatialAIM against\nthe annotation"
HEAD_MAPS = {"annotation": "Kidney, mapped\nannotation (K = {k})",
             "spatialaim": "Kidney, mapped\nSpatialAIM (K = {k})"}
# Under each map: its neighbourhood enrichment, then its key, both smaller than
# a caption.
FS_KEY = 32
KEY_DOT = 22
KEY_INDENT = 32       # dot to label
KEY_GAP = 26          # between two entries on a line
KEY_PITCH = FS_KEY * 1.25
MAP_GAP = 10          # map bottom to the metric line
ZOOM_GAP = 24         # zoom to the key beside it

ROWS = [
    {
        "sketch": ("sketch_spatial",
                   "Spatial organisation of\nthe mapped cell types"),
        "panels": [("plaus_nhood_mean_self_zscore",
                    "Mean z-score of\nneighbourhood enrichment diagonal"),
                   ("plaus_local_purity_zscore",
                    "Local spatial purity z-score")],
    },
    {
        "sketch": ("sketch_coherence",
                   "Transcriptional coherence of\nthe mapped cell types"),
        "panels": [("plaus_modularity_st_expression",
                    "Modularity over mapped\ncell types in spatial data"),
                   ("plaus_modularity_shared",
                    "Modularity over\nthe reference's cells")],
    },
]

# set by main(), once --mapper says which figure this is
drawing: Drawing


def text(body: str, x, y, size=FS_HEAD, colour=None, centre=False):
    return drawing.text(body, x, y, size=size,
                        align="centre" if centre else "left",
                        **({} if colour is None else {"colour": colour}))


def image(name: str, x, y, manifest: dict, folder=None) -> tuple[float, float]:
    element = drawing.image(name, x, y, manifest, folder)
    return element["width"], element["height"]


def key_block(x: float, y: float, width: float, entries: list[dict]) -> float:
    """A map's key: a dot and a label per entry, flowed left to right and
    wrapped onto a new line when ``width`` is used up; a label wider than
    ``width`` on its own is broken over lines. Returns the y below the key."""
    cursor, line, line_h = x, y, KEY_PITCH
    for entry in entries:
        label = wrap(entry["label"], width - KEY_INDENT, FS_KEY)
        w = KEY_INDENT + text_size(label, FS_KEY)[0]
        if cursor > x and cursor + w > x + width:
            cursor, line, line_h = x, line + line_h, KEY_PITCH
        drawing.dot(cursor + KEY_DOT / 2, line + FS_KEY * 0.6, KEY_DOT, entry["colour"])
        text(label, cursor + KEY_INDENT, line, size=FS_KEY)
        line_h = max(line_h, text_size(label, FS_KEY)[1] + FS_KEY * 0.25)
        cursor += w + KEY_GAP
    return line + line_h


def top_row(manifest: dict, kidney: dict, letters, y: float) -> float:
    """a, b: the kidney mapped as annotated and by SpatialAIM, each with its
    neighbourhood enrichment underneath, then its zoom (the boxed window) and its
    key side by side; c: the win counts. The three
    letters and captions share one line, the images the line below the tallest
    caption, and the slack is split into two equal gaps. Returns the y of the
    next row's letter line."""
    folder = kidney["folder"]
    w_wins = manifest["panels"]["wins"]["w"]
    w_map = kidney["panels"]["map_annotation"]["w"]
    gap = (WIDTH - 2 * PAD - w_wins - 2 * w_map) / 2
    xs = [PAD, PAD + w_map + gap, PAD + 2 * w_map + 2 * gap]
    sides = ["annotation", "spatialaim"]
    heads = [HEAD_MAPS[s].format(k=kidney[s]["k"]) for s in sides] + [HEAD_WINS]
    top = y + max(h.count("\n") + 1 for h in heads) * LINE + HEAD_GAP
    for head, x in zip(heads, xs):
        text(next(letters), x, y, size=FS_LETTER)
        text(head, x + HEAD_DX, y + 9)

    _w, h = image("wins", xs[2], top, manifest)
    bottom = top + h
    for side, x in zip(sides, xs[:2]):
        _w, h_map = image(f"map_{side}", x, top, kidney, folder)
        z = kidney[side]["metrics"]["nhood_mean_self_zscore"]
        line = top + h_map + MAP_GAP
        text(f"Neighbourhood enrichment z = {z:.0f}", x, line, size=FS_KEY, colour=GREY)
        below = line + KEY_PITCH * 1.3
        # the zoom (the window boxed on the map) at the left, the key beside it
        if kidney.get("zoom"):
            w_zoom, h_zoom = image(f"zoom_{side}", x, below, kidney, folder)
            key_x, key_w = x + w_zoom + ZOOM_GAP, w_map - w_zoom - ZOOM_GAP
            bottom = max(bottom, below + h_zoom)
        else:
            key_x, key_w = x, w_map
        bottom = max(bottom, key_block(key_x, below, key_w, kidney[side]["key"]))
    return bottom + ROW_GAP


# --------------------------------------------------------------------------
# the figure
# --------------------------------------------------------------------------

def build(manifest: dict, kidney: dict | None) -> None:
    panels = manifest["panels"]
    letters = iter("abcdefghijklmnopqrstuvwxyz")
    w_sketch = panels["sketch_spatial"]["w"]
    w_scatter = panels["plaus_nhood_mean_self_zscore"]["w"]
    h_row = panels["plaus_nhood_mean_self_zscore"]["h"]

    # Three cells per row on one baseline; the slack is split into two equal
    # gaps, so the row fills the figure's width exactly.
    gap = (WIDTH - 2 * PAD - w_sketch - 2 * w_scatter) / 2
    xs = [PAD, PAD + w_sketch + gap, PAD + w_sketch + w_scatter + 2 * gap]
    widths = [w_sketch, w_scatter, w_scatter]

    # ---------------- the top row (optional) ----------------
    y = ROW_TOP
    if kidney is not None:
        y = top_row(manifest, kidney, letters, y)

    for row in ROWS:
        cells = [row["sketch"], *row["panels"]]
        # each caption wrapped to the width its own cell leaves beside the letter
        heads = [wrap(body, width - HEAD_DX, FS_HEAD)
                 for (_name, body), width in zip(cells, widths)]
        # Every caption in the row *starts* on the same line, so the three
        # letters sit on one horizontal row; a caption shorter than the tallest
        # one leaves the slack above the image rather than above itself.
        lines = max(head.count("\n") + 1 for head in heads)
        top = y + lines * LINE + HEAD_GAP

        for head, x in zip(heads, xs):
            text(next(letters), x, y, size=FS_LETTER)
            text(head, x + HEAD_DX, y + 9)

        _w_sketch, h_sketch = image(row["sketch"][0], xs[0], top, manifest)
        # The two words naming the sketch's halves, on the empty strip each rect
        # leaves above it. Scaled by the sketch's own height -- it is shorter
        # than a scatter, so the row height would drift both labels downwards.
        for label, fraction in zip(("Observed", "Chance"),
                                   manifest["sketch_tops"][row["sketch"][0]]):
            text(label, xs[0] + 6,
                 top + fraction * h_sketch - text_size(label, FS_HEAD)[1])

        for x, (name, _body) in zip(xs[1:], row["panels"]):
            image(name, x, top, manifest)

        y = top + h_row + ROW_GAP

    # one colour bar for all four scatters, centred under the block
    bar = panels["dk_bar"]
    image("dk_bar", (WIDTH - bar["w"]) / 2, y - ROW_GAP + BAR_GAP, manifest)
    bottom = y - ROW_GAP + BAR_GAP + bar["h"] + PAD
    logger.info("canvas %.0f x %.0f units  (%.2f x %.2f in)",
                WIDTH, bottom, WIDTH / SCALE, bottom / SCALE)


def main() -> None:
    global drawing
    ap = argparse.ArgumentParser(description=__doc__)
    which = ap.add_mutually_exclusive_group(required=True)
    which.add_argument("--mapper", choices=MAPPERS,
                       help="whose panels to assemble, from figure_4/<mapper>/panels/")
    which.add_argument("--all-mappers", action="store_true",
                       help="the version aggregated over all mappers, figure_4/all_mappers/")
    ap.add_argument("--top-row", action=argparse.BooleanOptionalAction, default=None,
                    help="include the top row, win counts and kidney example (default: "
                         "only for " + ", ".join(sorted(TOP_ROW_DEFAULT)) + ")")
    args = ap.parse_args()
    mapper = args.mapper or "all_mappers"
    with_top = mapper in TOP_ROW_DEFAULT if args.top_row is None else args.top_row
    panels = paths.panel_dir(FIGURE, mapper)
    drawing = Drawing("f24", panels, paths.VAULT)

    manifest = json.loads((panels / "panels.json").read_text(encoding="utf-8"))
    kidney = None
    if with_top:
        path = paths.figure_dir(FIGURE, mapper) / "kidney_example" / "panels.json"
        if not path.exists():
            raise SystemExit(f"{path} missing - run 04_kidney_example/ first (its README)")
        kidney = {**json.loads(path.read_text(encoding="utf-8")), "folder": path.parent}
    build(manifest, kidney)
    drawing.check()
    drawing.write(paths.drawing(FIGURE, mapper))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
