"""Building an Excalidraw document, which is how every figure of the paper is
assembled.

The division of labour is the same in all four figures: ``panels.py`` writes the
data as PNGs, and the drawing adds what makes them readable -- panel letters,
headings, arrows and legends. The drawing embeds each panel by wikilink, so
re-running ``panels.py`` updates the figure in place, and the document stays
editable by hand in Obsidian afterwards.

This module is the part every ``make_excalidraw.py`` had copied: element
construction, text measurement, and writing the plugin's file format. What each
figure does with them -- its rows, its gaps, its headings -- stays in its own
section.

Sizes are in canvas units, ``SCALE`` per print inch. The paper places a figure
at its 5.17 in text width, so a 7.2 in drawing is scaled by 0.71 and every point
size is authored 1/0.71 larger than it should print.
"""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path

from .style import GREY, INK, SCALE

logger = logging.getLogger(__name__)

FONT = 2              # Helvetica, as in figure 1
WIDTH = 7.2 * SCALE   # the width every figure is drawn at

# Helvetica is proportional, so a text element's width has to be estimated per
# character. One average glyph is not enough: it under-measures capital-heavy
# words -- "WANN" and "TACCO" ran into their neighbours in the key -- and
# over-measures narrow ones. These are rough advance widths in em.
NARROW, WIDE = "iljtfI.,;:!|'\"()[]{}- ", "MWmw"


def advance(ch: str) -> float:
    if ch in NARROW:
        return 0.30
    if ch in WIDE:
        return 0.87
    return 0.70 if ch.isupper() or ch.isdigit() else 0.54


def text_size(body: str, size: float) -> tuple[float, float]:
    lines = body.split("\n")
    widest = max(sum(advance(c) for c in line) for line in lines)
    return widest * size, 1.2 * size * len(lines)


def wrap(body: str, width: float, size: float) -> str:
    """Break text at spaces so that no line is wider than ``width``. Line breaks
    already in the text are kept: those are where someone wanted them."""
    lines = []
    for paragraph in body.split("\n"):
        line = ""
        for word in paragraph.split():
            candidate = f"{line} {word}".strip()
            if line and text_size(candidate, size)[0] > width:
                lines.append(line)
                line = word
            else:
                line = candidate
        lines.append(line)
    return "\n".join(lines)


class Drawing:
    """One Excalidraw document under construction.

    ``prefix`` goes into every element id, so the ids of two figures never
    collide if a reader copies elements between them. ``panel_dir`` is the folder
    holding the PNGs; it must lie inside ``vault``, because the plugin resolves an
    embedded file by its vault-relative wikilink.
    """

    def __init__(self, prefix: str, panel_dir: Path, vault: Path,
                 width: float = WIDTH, font: int = FONT):
        self.prefix = prefix
        self.panel_dir = Path(panel_dir)
        self.vault = Path(vault)
        self.width = width
        self.font = font
        self.elements: list[dict] = []
        self.embedded: dict[str, str] = {}

    # -- elements ----------------------------------------------------------
    def element(self, kind: str, x: float, y: float, w: float, h: float,
                **kw) -> dict:
        element = {
            "id": f"{self.prefix}{kind[:4]}{len(self.elements):03d}",
            "type": kind,
            "x": round(x, 2), "y": round(y, 2),
            "width": round(w, 2), "height": round(h, 2),
            "angle": 0,
            "strokeColor": kw.pop("strokeColor", INK),
            "backgroundColor": kw.pop("backgroundColor", "transparent"),
            "fillStyle": "solid",
            "strokeWidth": kw.pop("strokeWidth", 1),
            "strokeStyle": "solid",
            "roughness": 0,
            "opacity": 100,
            "groupIds": [],
            "frameId": None,
            "roundness": kw.pop("roundness", None),
            "seed": 1000 + len(self.elements),
            "version": 1,
            "versionNonce": 2000 + len(self.elements),
            "isDeleted": False,
            "boundElements": [],
            "updated": 1755600000000,
            "link": None,
            "locked": False,
            "index": f"a{len(self.elements):04d}",
        }
        element.update(kw)
        self.elements.append(element)
        return element

    def text(self, body: str, x: float, y: float, size: float = 40,
             colour: str = INK, align: str = "left",
             box: float | None = None) -> dict:
        """``align`` says what ``x`` means: the left edge, the right edge, or the
        centre of the text block.

        Without ``box`` the element is as wide as ``text_size`` estimates, and
        that estimate is what gets aligned -- good enough for a single label, but
        a column of them comes out ragged, because the estimate is off by a
        different amount per word. ``box`` fixes the element's width instead and
        lets Excalidraw align the glyphs inside it, which is flush. It has to be
        at least as wide as the widest label, or that label wraps.
        """
        w, h = text_size(body, size)
        width = w if box is None else box
        left = {"left": x, "right": x - width, "centre": x - width / 2}[align]
        inside = {"left": "left", "right": "right", "centre": "center"}[align]
        return self.element("text", left, y, width, h, strokeColor=colour,
                            fontSize=size, fontFamily=self.font, text=body,
                            originalText=body, rawText=body,
                            textAlign="left" if box is None else inside,
                            verticalAlign="top", containerId=None,
                            autoResize=box is None, lineHeight=1.2)

    def image(self, name: str, x: float, y: float, manifest: dict,
              folder: Path | None = None) -> dict:
        """A panel PNG, embedded by wikilink so re-exporting updates the
        drawing. Placed at the size ``panels.json`` recorded, not at its pixel
        count. ``folder`` overrides ``panel_dir`` for an image kept elsewhere in
        the vault, such as a screenshot taken by hand."""
        panel = manifest["panels"][name]
        file_id = hashlib.sha1(name.encode()).hexdigest()
        folder = self.panel_dir if folder is None else Path(folder)
        link = f"{folder.relative_to(self.vault).as_posix()}/{panel['file']}"
        self.embedded[file_id] = link
        return self.element("image", x, y, panel["w"], panel["h"],
                            fileId=file_id, status="saved", scale=[1, 1],
                            crop=None)

    def dot(self, cx: float, cy: float, diameter: float, colour: str) -> dict:
        return self.element("ellipse", cx - diameter / 2, cy - diameter / 2,
                            diameter, diameter, strokeColor=colour,
                            backgroundColor=colour)

    def bar(self, cx: float, cy: float, w: float, h: float, colour: str) -> dict:
        """A rounded vertical bar, centred on (cx, cy) -- the mapper swatch."""
        return self.element("rectangle", cx - w / 2, cy - h / 2, w, h,
                            strokeColor=colour, backgroundColor=colour,
                            roundness={"type": 3})

    def line(self, x: float, y: float, length: float, colour: str,
             stroke: float = 1) -> dict:
        """A horizontal rule, the legend swatch for a plotted curve."""
        return self.element("line", x, y, length, 0, strokeColor=colour,
                            strokeWidth=stroke, points=[[0, 0], [length, 0]],
                            lastCommittedPoint=None, startBinding=None,
                            endBinding=None, startArrowhead=None,
                            endArrowhead=None, polygon=False)

    def arrow(self, points: list[tuple[float, float]], colour: str = GREY,
              **kw) -> dict:
        """Points are absolute; corners are rounded, as in figure 1."""
        x0, y0 = points[0]
        relative = [[round(x - x0, 2), round(y - y0, 2)] for x, y in points]
        xs = [p[0] for p in relative]
        ys = [p[1] for p in relative]
        return self.element("arrow", x0, y0, max(xs) - min(xs), max(ys) - min(ys),
                            strokeColor=colour, roundness={"type": 2},
                            points=relative, lastCommittedPoint=None,
                            startBinding=None, endBinding=None,
                            startArrowhead=None, endArrowhead="arrow", **kw)

    # -- output ------------------------------------------------------------
    def check(self) -> None:
        """Warn about anything that leaves the figure, which a longer heading
        does easily -- the drawing itself will not complain."""
        for element in self.elements:
            right = element["x"] + element["width"]
            if element["x"] < 0 or right > self.width:
                logger.warning("%s sticks out by %.0f units: %r", element["type"],
                               max(-element["x"], right - self.width),
                               element.get("text", element["id"]))

    def write(self, path: Path) -> None:
        drawing = {
            "type": "excalidraw",
            "version": 2,
            "source": "https://github.com/zsviczian/obsidian-excalidraw-plugin",
            "elements": self.elements,
            "appState": {"theme": "light", "viewBackgroundColor": "#ffffff",
                         "currentItemFontFamily": self.font,
                         "currentItemFontSize": 40,
                         "currentItemRoughness": 0, "gridSize": None},
            "files": {},
        }
        body = [
            "---", "", "excalidraw-plugin: parsed", "tags: [excalidraw]", "", "---",
            "==⚠  Switch to EXCALIDRAW VIEW in the MORE OPTIONS menu of this"
            " document. ⚠==", "", "", "# Excalidraw Data", "",
            "## Text Elements",
        ]
        for element in self.elements:
            if element["type"] == "text":
                body.append(f"{element['text']} ^{element['id']}\n")
        body += ["", "## Embedded Files"]
        for file_id, link in self.embedded.items():
            body.append(f"{file_id}: [[{link}]]\n")
        body += ["", "%%", "## Drawing", "```json", json.dumps(drawing), "```", "%%"]
        Path(path).write_text("\n".join(body), encoding="utf-8")
        logger.info("%d elements, %d embedded panels -> %s", len(self.elements),
                    len(self.embedded), path)
