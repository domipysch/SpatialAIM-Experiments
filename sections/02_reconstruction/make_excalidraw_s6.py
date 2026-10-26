# Assemble supplementary figure 6, the peak memory of the K sweep from the Leiden
# start (moved out of figure 2), as an Excalidraw drawing in the Obsidian vault.
#
# One panel, written by panels.py into figure_S6/panels/ and embedded by wikilink,
# with its heading above and the key to the five mapper colours stacked to its
# right. Sizes as in figure 2 (authored 1/0.71 larger than they print).

import json
import logging

from spatialaim_experiments import paths
from spatialaim_experiments.excalidraw import Drawing
from spatialaim_experiments.style import MAPPER_COLOURS, MAPPER_LABELS, MAPPERS

logger = logging.getLogger(__name__)

FIGURE = "S6"
PANELS = paths.panel_dir(FIGURE)
OUT = paths.drawing(FIGURE)

PAD = 14
FS_HEAD = 40
FS_LEG = 40
MARK_W, MARK_H = 16, 48   # the mapper swatch, one line high, as in figure 2's key
LEG_INDENT = MARK_W + 18
LINE = FS_LEG * 1.2
HEADING = "Peak memory, K sweep from the Leiden start"

drawing = Drawing("fs6", PANELS, paths.VAULT)


def build(manifest: dict) -> None:
    panel = manifest["panels"]["memory"]
    drawing.text(HEADING, PAD, PAD, size=FS_HEAD)
    top = PAD + FS_HEAD * 1.2 + 20
    drawing.image("memory", PAD, top, manifest)
    # the key, one mapper per line, level with the plotting box
    x = PAD + panel["w"] + 40
    for i, mapper in enumerate(MAPPERS):
        y = top + 30 + i * (LINE + 8)
        drawing.bar(x + MARK_W / 2, y + FS_LEG * 0.6, MARK_W, MARK_H, MAPPER_COLOURS[mapper])
        drawing.text(MAPPER_LABELS[mapper], x + LEG_INDENT, y, size=FS_LEG)


def main() -> None:
    manifest = json.loads((PANELS / "panels.json").read_text(encoding="utf-8"))
    build(manifest)
    drawing.check()
    drawing.write(OUT)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
