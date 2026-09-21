"""Render an Excalidraw results figure to the PDF the paper includes.

    python -m aim_experiments.render 2_3
        00_Obsidian/Excalidraw/figure_2_3.excalidraw.md
            -> AIM-Overleaf/figures/fig_2_3.pdf
    python -m aim_experiments.render 2_3 --out preview.png   (a look, not the paper)

It reads the drawing itself, not the script that generates it, so anything moved
or edited by hand in Obsidian is what ends up in the paper. This is a stand-in
for the plugin's own "export image": the panels are embedded at their stored
resolution and everything else is drawn as vector, but the text is set in Arial
rather than Excalidraw's Helvetica -- close in metrics, not identical. For the
camera-ready figure, export from Excalidraw itself.
"""

import argparse
import json
import logging
import re
from pathlib import Path

import matplotlib as mpl
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.image import imread
from matplotlib.patches import Ellipse, FancyArrow, FancyBboxPatch, Rectangle

from . import paths
from .style import SCALE

logger = logging.getLogger(__name__)
logging.getLogger("fontTools").setLevel(logging.WARNING)  # it narrates subsetting

VAULT = paths.VAULT

MARGIN = 8       # trimmed like Excalidraw's own export, with a little air
FONT = "Arial"   # what Windows ships in place of Excalidraw's Helvetica

# Excalidraw's three stroke styles, in matplotlib dash units
DASHES = {"solid": "solid", "dashed": (0, (5.5, 4.0)), "dotted": (0, (1.0, 3.0))}

mpl.rcParams.update({"pdf.fonttype": 42, "font.family": "sans-serif",
                     "font.sans-serif": [FONT, "Helvetica", "DejaVu Sans"]})

# --------------------------------------------------------------------------
# lz-string decompressFromBase64: the plugin rewrites the drawing block in that
# form whenever the file is saved from the Excalidraw view
# --------------------------------------------------------------------------
KEY = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/="


def decompress(compressed: str) -> str:
    data = "".join(compressed.split())
    index = {c: i for i, c in enumerate(KEY)}
    dictionary, out = [0, 1, 2], []
    enlarge, size, bits_per = 4, 4, 3
    value, position, pointer = index[data[0]], 32, 1

    def read(count: int) -> int:
        nonlocal value, position, pointer
        result, maximum, power = 0, 1 << count, 1
        while power != maximum:
            bit = value & position
            position >>= 1
            if position == 0:
                position = 32
                value = index[data[pointer]]
                pointer += 1
            result |= (1 if bit > 0 else 0) * power
            power <<= 1
        return result

    marker = read(2)
    character = chr(read(8 if marker == 0 else 16))
    dictionary.append(character)
    word = character
    out.append(character)
    while True:
        if pointer > len(data):
            return ""
        marker = read(bits_per)
        if marker in (0, 1):
            dictionary.append(chr(read(8 if marker == 0 else 16)))
            marker, size, enlarge = size, size + 1, enlarge - 1
        elif marker == 2:
            return "".join(out)
        if enlarge == 0:
            enlarge, bits_per = 1 << bits_per, bits_per + 1
        entry = dictionary[marker] if marker < len(dictionary) else word + word[0]
        out.append(entry)
        dictionary.append(word + entry[0])
        size, enlarge = size + 1, enlarge - 1
        word = entry
        if enlarge == 0:
            enlarge, bits_per = 1 << bits_per, bits_per + 1


def load(drawing_file: Path) -> tuple[list[dict], dict[str, Path]]:
    text = drawing_file.read_text(encoding="utf-8")
    packed = re.search(r"```compressed-json\n(.*?)\n```", text, re.S)
    plain = re.search(r"```json\n(.*?)\n```", text, re.S)
    drawing = json.loads(decompress(packed.group(1)) if packed else plain.group(1))
    # the plugin shortens the wikilinks to bare file names when it saves, so a
    # link is resolved the way Obsidian resolves it: by name, anywhere below the
    # vault
    files = {}
    for file_id, link in re.findall(r"^([0-9a-f]{40}): \[\[(.+?)\]\]", text, re.M):
        path = VAULT / link
        if not path.exists():
            path = next(VAULT.rglob(Path(link).name))
        files[file_id] = path
    logger.info("%s: %d elements, %d panels", "compressed" if packed else "plain",
                len(drawing["elements"]), len(files))
    return drawing["elements"], files


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("figure", nargs="?", default="2_2",
                    help="which drawing to render, e.g. 2_3 or 2_2_wann "
                         "(default 2_2)")
    ap.add_argument("--out", type=Path,
                    help="write here instead of AIM-Overleaf/figures/fig_<figure>.pdf; "
                         "give it a .png suffix for a quick look that never touches "
                         "the paper")
    args = ap.parse_args()
    out = args.out or paths.PAPER_FIGURES / f"fig_{args.figure}.pdf"

    elements, files = load(paths.EXCALIDRAW / f"figure_{args.figure}.excalidraw.md")
    x0 = min(e["x"] for e in elements) - MARGIN
    y0 = min(e["y"] for e in elements) - MARGIN
    x1 = max(e["x"] + e["width"] for e in elements) + MARGIN
    y1 = max(e["y"] + e["height"] for e in elements) + MARGIN

    fig = plt.figure(figsize=((x1 - x0) / SCALE, (y1 - y0) / SCALE))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)   # Excalidraw's y grows downwards
    ax.axis("off")

    for element in sorted(elements, key=lambda e: e.get("index", "")):
        kind = element["type"]
        x, y = element["x"], element["y"]
        w, h = element["width"], element["height"]
        colour = element["strokeColor"]
        fill = element["backgroundColor"]
        if kind == "image":
            ax.imshow(imread(files[element["fileId"]]),
                      extent=[x, x + w, y + h, y], zorder=2, interpolation="none")
        elif kind == "rectangle":
            face = "none" if fill == "transparent" else fill
            if element.get("roundness"):
                # Excalidraw's adaptive radius, near enough: a quarter of the
                # shorter side. Without this a rounded swatch renders square and
                # the preview misreports the drawing.
                radius = min(w, h) * 0.25
                ax.add_patch(FancyBboxPatch(
                    (x + radius, y + radius), w - 2 * radius, h - 2 * radius,
                    boxstyle=f"round,pad={radius},rounding_size={radius}",
                    ec=colour, lw=element["strokeWidth"], fc=face))
            else:
                ax.add_patch(Rectangle((x, y), w, h, ec=colour,
                                       lw=element["strokeWidth"], fc=face))
        elif kind == "ellipse":
            ax.add_patch(Ellipse((x + w / 2, y + h / 2), w, h, ec=colour,
                                 lw=element["strokeWidth"] * 0.8,
                                 fc="none" if fill == "transparent" else fill))
        elif kind in ("line", "arrow"):
            points = np.asarray(element["points"], dtype=float) + (x, y)
            # A closed line with a background is a filled polygon in Excalidraw
            # -- that is how a shape it has no primitive for, such as a star, is
            # drawn.
            if fill != "transparent" and element.get("polygon"):
                ax.fill(points[:, 0], points[:, 1], facecolor=fill,
                        edgecolor=colour, lw=element["strokeWidth"], zorder=3)
            ax.plot(points[:, 0], points[:, 1], color=colour, lw=element["strokeWidth"],
                    linestyle=DASHES.get(element.get("strokeStyle"), "solid"),
                    solid_capstyle="round", zorder=3)
            if element.get("endArrowhead") == "arrow":
                start, end = points[-2], points[-1]
                step = end - start
                length = float(np.hypot(*step)) or 1.0
                head = 13.0
                ax.add_patch(FancyArrow(*(end - step / length * head), *(step / length * head),
                                        width=0, head_width=9, head_length=head,
                                        color=colour, length_includes_head=True, zorder=3))
        elif kind == "text":
            size = element["fontSize"] * 72 / SCALE
            angle = float(element.get("angle") or 0.0)
            if angle:
                # Excalidraw turns an element about its own centre, and its angle
                # runs clockwise in radians -- matplotlib's runs anticlockwise in
                # degrees. Rotated labels (the "scores" axis title) came out
                # horizontal and looked as if they overlapped the data.
                ax.text(x + w / 2, y + h / 2, element["text"], color=colour,
                        va="center", ha="center", rotation=-np.degrees(angle),
                        rotation_mode="anchor", fontsize=size, linespacing=1.2,
                        zorder=4)
            else:
                ax.text(x, y, element["text"], color=colour, va="top", ha="left",
                        fontsize=size, linespacing=1.2, zorder=4)

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=SCALE)
    logger.info("%.2f x %.2f in -> %s", (x1 - x0) / SCALE, (y1 - y0) / SCALE, out)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
