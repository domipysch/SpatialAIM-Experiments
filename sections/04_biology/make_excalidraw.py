# Assemble the Results 2.4 figure -- biological and topological plausibility --
# as an Excalidraw drawing in the Obsidian vault.
#
# Same construction as figure_2_2 and figure_2_3: the data panels are the PNGs
# written by panels.py and embedded by wikilink, so re-running that script
# updates them in place. This script draws only what makes them readable --
# letters, headings and the two words naming each sketch half. Anything else
# belongs in the caption.
#
# One mapper per invocation (``--mapper``, default nearest_centroid, the figure
# the paper prints); every other mapper is assembled from figure_2_4_<mapper>/
# into figure_2_4_<mapper>.excalidraw.md.
#
# Two rows of three, one row per criterion, and every cell its own panel:
#
#   a sketch | b neighbourhood enrichment | c local spatial purity
#   d sketch | e modularity over the slice | f modularity over the reference
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
from spatialaim_experiments.style import MAPPERS, SCALE

logger = logging.getLogger(__name__)

FIGURE = "2_4"
# the mapper panels.py draws the figure itself with
FIGURE_MAPPER = "nearest_centroid"
# set by main() from --mapper; the figure's own mapper keeps the unsuffixed
# names, every other one gets its copy of the same drawing beside it
PANELS = paths.panel_dir(FIGURE)
OUT = paths.drawing(FIGURE)

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
ROWS = [
    {
        "sketch": ("a", "sketch_spatial",
                   "Spatial organisation of\nthe mapped cell types"),
        "panels": [("b", "plaus_nhood_mean_self_zscore",
                    "Mean z-score of\nneighbourhood enrichment diagonal"),
                   ("c", "plaus_local_purity_zscore",
                    "Local spatial purity z-score")],
    },
    {
        "sketch": ("d", "sketch_coherence",
                   "Transcriptional coherence of\nthe mapped cell types"),
        "panels": [("e", "plaus_modularity_st_expression",
                    "Modularity over mapped\ncell types in spatial data"),
                   ("f", "plaus_modularity_shared",
                    "Modularity over\nthe reference's cells")],
    },
]

# set by main(), once --mapper says which figure this is
drawing = Drawing("f24", PANELS, paths.VAULT)


def text(body: str, x, y, size=FS_HEAD, colour=None, centre=False):
    return drawing.text(body, x, y, size=size,
                        align="centre" if centre else "left",
                        **({} if colour is None else {"colour": colour}))


def image(name: str, x, y, manifest: dict) -> tuple[float, float]:
    element = drawing.image(name, x, y, manifest)
    return element["width"], element["height"]


# --------------------------------------------------------------------------
# the figure
# --------------------------------------------------------------------------

def build(manifest: dict) -> None:
    panels = manifest["panels"]
    w_sketch = panels["sketch_spatial"]["w"]
    w_scatter = panels["plaus_nhood_mean_self_zscore"]["w"]
    h_row = panels["plaus_nhood_mean_self_zscore"]["h"]

    # Three cells per row on one baseline; the slack is split into two equal
    # gaps, so the row fills the figure's width exactly.
    gap = (WIDTH - 2 * PAD - w_sketch - 2 * w_scatter) / 2
    xs = [PAD, PAD + w_sketch + gap, PAD + w_sketch + w_scatter + 2 * gap]
    widths = [w_sketch, w_scatter, w_scatter]

    y = ROW_TOP
    for row in ROWS:
        cells = [row["sketch"], *row["panels"]]
        # each caption wrapped to the width its own cell leaves beside the letter
        heads = [wrap(body, width - HEAD_DX, FS_HEAD)
                 for (_letter, _name, body), width in zip(cells, widths)]
        # Every caption in the row *starts* on the same line, so the three
        # letters sit on one horizontal row; a caption shorter than the tallest
        # one leaves the slack above the image rather than above itself.
        lines = max(head.count("\n") + 1 for head in heads)
        top = y + lines * LINE + HEAD_GAP

        for (letter, _name, _body), head, x in zip(cells, heads, xs):
            text(letter, x, y, size=FS_LETTER)
            text(head, x + HEAD_DX, y + 9)

        _w_sketch, h_sketch = image(row["sketch"][1], xs[0], top, manifest)
        # The two words naming the sketch's halves, on the empty strip each rect
        # leaves above it. Scaled by the sketch's own height -- it is shorter
        # than a scatter, so the row height would drift both labels downwards.
        for label, fraction in zip(manifest["sketch_labels"],
                                   manifest["sketch_tops"][row["sketch"][1]]):
            text(label, xs[0] + 6,
                 top + fraction * h_sketch - text_size(label, FS_HEAD)[1])

        for x, (_letter, name, _body) in zip(xs[1:], row["panels"]):
            image(name, x, top, manifest)

        y = top + h_row + ROW_GAP

    # one colour bar for all four scatters, centred under the block
    bar = panels["dk_bar"]
    image("dk_bar", (WIDTH - bar["w"]) / 2, y - ROW_GAP + BAR_GAP, manifest)
    bottom = y - ROW_GAP + BAR_GAP + bar["h"] + PAD
    logger.info("canvas %.0f x %.0f units  (%.2f x %.2f in)",
                WIDTH, bottom, WIDTH / SCALE, bottom / SCALE)


def main() -> None:
    global PANELS, OUT, drawing
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mapper", default=FIGURE_MAPPER, choices=MAPPERS,
                    help=f"whose panels to assemble (default {FIGURE_MAPPER}, the "
                         "figure itself)")
    mapper = ap.parse_args().mapper
    copy = None if mapper == FIGURE_MAPPER else mapper
    PANELS = paths.panel_dir(FIGURE, copy)
    OUT = paths.drawing(FIGURE, copy)
    drawing = Drawing("f24", PANELS, paths.VAULT)

    manifest = json.loads((PANELS / "panels.json").read_text(encoding="utf-8"))
    build(manifest)
    drawing.check()
    drawing.write(OUT)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
