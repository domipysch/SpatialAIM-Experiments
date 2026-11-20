# Assemble the Results 2.5 GUI use-case figure (Oliveira et al. 2025 CRC P2) as an
# Excalidraw drawing in the Obsidian vault:
#
#   a  the GUI's Explore page at K = 10, numbered callouts
#   b  Oliveira et al. vs SpatialAIM, two vertical lanes
#   c  the cell-type modal's reference-UMAP maps of SPP1 and STAB1
#   d  cell area of cell types 2 vs 3 (panels.py, from the GUI's shape-features card)
#   e  K = 10, K = 16 and the static Level1 annotation as identical pairs
#      (K slider, reference UMAP, Xenium mapping) in one row
#
# Images are GUI screenshots from screenshots.py (figure_5/screenshots/) plus the
# area plot, both from panels.py (figure_5/panels/); the numbers quoted come from
# figure_5/panels/panels.json. Letters, headings, callouts and the schematic are
# drawn here.
#
#   python make_excalidraw.py            -> figure_5.excalidraw.md

import json
import logging
from string import ascii_lowercase

from spatialaim_experiments import paths
from spatialaim_experiments.excalidraw import WIDTH, Drawing, text_size
from spatialaim_experiments.style import INK, png_size

logger = logging.getLogger(__name__)

FIGURE = "5"
PANELS = paths.panel_dir(FIGURE)
OUT = paths.drawing(FIGURE)
SHOTS = PANELS.parent / "screenshots"
AREA_DIR = PANELS   # panels.py's cell-area plot

PAD = 14
GAP = 40
FS_LETTER = 50
FS_HEAD = 40
FS_SMALL = 34
FS_BOX = 32
FS_NOTE = 30
HEAD_DX = 44
LETTER_DROP = FS_LETTER * 1.2 + 6
ROW_GAP = 40
SPAN = WIDTH - 2 * PAD

PAPER_C, AIM_C = "#868e96", "#0072B2"
RESULT_FILL = "#fff4e6"
HINT_C = "#d6336c"   # interactivity hints, away from the GUI's state colours

# (a)'s callouts, explained in the caption (positions: screenshots.py CALLOUT_TARGETS):
# 1 select data, 2 configure method / run, 3 starting clusters + mapper,
# 4 K slider (dots: Pareto-optimal K), 5 previous / next Pareto K, best K,
# 6 Choose K card (criteria per K; open in the supplement),
# 7 cell types (click: hide, right click: inspect), 8 reference UMAP,
# 9 mapped Xenium cells

# hand placements (2026-11-19, in Obsidian), drawing units from the computed spot:
# 1 and 2 off the sidebar's left edge, 5 beside Best K, 6 mid Choose K card
CALLOUT_NUDGE = {1: (37.6, -3.3), 2: (26.2, -18.0), 5: (88.4, 91.6), 6: (634.3, -7.0)}

drawing: Drawing = None
letters = iter(ascii_lowercase)


# --------------------------------------------------------------------------
# primitives
# --------------------------------------------------------------------------

def text(body, x, y, size=FS_HEAD, **kw):
    return drawing.text(body, x, y, size=size, **kw)


def label(heading: str, x: float, top: float) -> None:
    """The next panel letter and its heading, on the line above ``top``."""
    text(next(letters), x, top - LETTER_DROP, size=FS_LETTER)
    text(heading, x + HEAD_DX, top - LETTER_DROP + 9)


def aspect(name: str, folder=None) -> float:
    w, h = png_size((folder or SHOTS) / f"{name}.png")
    return w / h


def shot(manifest: dict, name: str, x: float, y: float, width: float = None,
         height: float = None, folder=None, frame: bool = True) -> tuple[float, float]:
    """Place a screenshot (or, from ``folder``, a plot) at ``width`` or ``height``
    (aspect kept), framed."""
    folder = folder or SHOTS
    if width is None:
        width = height * aspect(name, folder)
    h = round(width / aspect(name, folder), 1)
    manifest["panels"][name] = {"file": f"{name}.png", "w": round(width, 1), "h": h}
    drawing.image(name, x, y, manifest, folder=folder)
    if frame:
        drawing.element("rectangle", x, y, width, h, strokeColor="#ced4da")
    return width, h


def box(x, y, w, h, body, colour, fill="#ffffff", size=FS_BOX):
    drawing.element("rectangle", x, y, w, h, strokeColor=colour, backgroundColor=fill,
                    strokeWidth=2, roundness={"type": 3})
    th = text_size(body, size)[1]
    text(body, x + w / 2, y + (h - th) / 2, size=size, align="centre", box=w - 16)


def callout(n: int, cx: float, cy: float, d: float = 62) -> None:
    drawing.dot(cx, cy, d, INK)
    size = d * 0.65
    th = text_size(str(n), size)[1]
    # the line box puts digits slightly high: nudge them to the dot's centre
    text(str(n), cx, cy - th / 2 + 0.07 * size, size=size, colour="#ffffff", align="centre", box=d)


def next_top(bottom: float) -> float:
    return bottom + ROW_GAP + LETTER_DROP


# --------------------------------------------------------------------------
# blocks: each draws itself (letter included) and returns its (width, height)
# below ``top``, the line its letter sits above
# --------------------------------------------------------------------------

def explore(m: dict, x: float, top: float, w: float):
    label("SpatialAIM GUI: Explore page", x, top)
    w, h = shot(m, "explore_k10", x, top, width=w)
    spec = json.loads((SHOTS / "callouts.json").read_text())
    cx, cy, cw, ch = spec["crop"]
    d = 62
    for n, (px, py) in spec["points"].items():
        # points sit 26 CSS px left of their target: keep the larger dot clear of it
        dx = max(x + (px - cx) / cw * w - (d / 2 - 20), x + d / 2 - PAD)   # may overhang into the margin
        ndx, ndy = CALLOUT_NUDGE.get(int(n), (0, 0))
        callout(int(n), dx + ndx, top + (py - cy) / ch * h + ndy, d=d)
    return w, h


ROUTES = [
    ("Oliveira et al.", PAPER_C, [(0, 1, "scRNA-seq\n+ Visium HD"), (1, 1, "label\nVisium HD bins"),
                                  (2, 1, "DGE on\nmacrophage bins"),
                                  (3, 1, "find substates\nSPP1+ / SELENOP+")]),
    ("SpatialAIM", AIM_C, [(0, 1, "scRNA-seq\n+ Xenium"), (1, 2, "run default\nSpatialAIM"),
                           (3, 1, "same substates")]),
]


def routes(m: dict, x: float, top: float, w: float, h: float):
    """The two routes, one column each; SpatialAIM's one run spans the steps the
    paper needed separately (its 50 um tumour-periphery step is left to the text).
    The boxes fill the height ``h``."""
    label("Finding Macrophage substates", x, top)
    head, gap, lane_gap = 50, 28, 30
    bh = (h - head - 3 * gap) / 4
    lane_w = (w - lane_gap) / 2
    for i, (title, colour, steps) in enumerate(ROUTES):
        lx = x + i * (lane_w + lane_gap)
        text(title, lx + lane_w / 2, top, size=FS_SMALL, colour=colour, align="centre", box=lane_w)
        prev = None
        for row, span, body in steps:
            by = top + head + row * (bh + gap)
            bhh = span * bh + (span - 1) * gap
            box(lx, by, lane_w, bhh, body, colour, fill=RESULT_FILL if row == 3 else "#ffffff",
                size=FS_BOX - 4)
            if prev is not None:
                drawing.arrow([(lx + lane_w / 2, prev + 4), (lx + lane_w / 2, by - 6)], colour=colour)
            prev = by + bhh
    return w, head + 4 * bh + 3 * gap


# the cell-type modal's gene maps, re-rendered by panels.py with larger axis text
MARKER_MAPS = (("genemap_spp1", "Cell state 2 · SPP1"), ("genemap_stab1", "Cell state 3 · STAB1"))


def markers(m: dict, x: float, top: float, h: float, sub: float = 40, gap: float = 20):
    """The SPP1 and STAB1 maps side by side at height ``h``, each under its name."""
    label("Macrophage substates: marker genes", x, top)
    cx = x
    for name, title in MARKER_MAPS:
        wm = (h - sub) * aspect(name, PANELS)
        text(title, cx + wm / 2, top, size=FS_NOTE, align="centre", box=wm)
        shot(m, name, cx, top + sub, height=h - sub, folder=PANELS, frame=False)
        cx += wm + gap
    return cx - gap - x, h


def area(m: dict, x: float, top: float, h: float, sub: float = 40):
    """Cell area of cell states 2 vs 3 at K = 10, its plot aligned with the marker maps."""
    label("Cell area", x, top)
    text("Cell states 2 vs 3", x, top, size=FS_NOTE)
    return shot(m, "area", x, top + sub, height=h - sub, folder=AREA_DIR, frame=False)


PAIRS = (("kslider_k10", "umap_k10", "spatial_k10", "K = 10 (best overall)"),
         ("kslider_k16", "umap_k16", "spatial_k16", "K = 16 (Pareto-optimal)"),
         (None, "umap_level1", "spatial_level1", "Level1 annotation (static)"))


def stepping(m: dict, x: float, top: float, w: float):
    """K = 10, K = 16 and the static Level1 annotation as identical pairs in one
    row: the K slider at that value (none for the annotation, which no K changes),
    the reference UMAP and the mapping onto the Xenium cells. Level1's two plots
    are the Leiden page's annotation row, coloured to match K = 16."""
    label("Stepping from K = 10 to the Pareto-optimal K = 16", x, top)
    pair_gap, inner, sub = 46, 8, 40
    w_pair = (w - 2 * pair_gap) / 3
    r_u, r_s, r_k = aspect("umap_k10"), aspect("spatial_k10"), aspect("kslider_k10")
    h_map = (w_pair - inner) / (r_u + r_s)
    h_k = w_pair / r_k
    ym = top + sub + h_k + 10
    for i, (slider, umap, spatial, title) in enumerate(PAIRS):
        px = x + i * (w_pair + pair_gap)
        text(title, px + w_pair / 2, top, size=FS_NOTE, align="centre", box=w_pair)
        if slider:
            shot(m, slider, px, top + sub, width=w_pair)
        shot(m, umap, px, ym, height=h_map)
        shot(m, spatial, px + h_map * r_u + inner, ym, height=h_map)
    # the step itself, from the K = 10 pair to the K = 16 pair
    x16, y_mid = x + w_pair + pair_gap, ym + h_map / 2
    drawing.arrow([(x16 - pair_gap + 6, y_mid), (x16 - 6, y_mid)], colour=HINT_C, strokeWidth=3)
    return w, ym + h_map - top


# --------------------------------------------------------------------------
# the figure
# --------------------------------------------------------------------------

def layout(m: dict) -> float:
    """GUI; routes beside the marker maps and the cell area; stepping. Returns
    the canvas height."""
    top = PAD + LETTER_DROP
    _, h = explore(m, PAD, top, SPAN)
    top = next_top(top + h)
    w_b, sub, gap = 600, 40, 20
    # marker maps and area at one height, filling what the routes leave; the
    # routes' boxes then fill that same height
    a_maps = sum(aspect(n, PANELS) for n, _ in MARKER_MAPS)
    a_area = aspect("area", AREA_DIR)
    hh = (SPAN - w_b - 2 * GAP - gap) / (a_maps + a_area) + sub
    _, hb = routes(m, PAD, top, w_b, hh)
    wm, hm = markers(m, PAD + w_b + GAP, top, h=hh, sub=sub, gap=gap)
    area(m, PAD + w_b + GAP + wm + GAP, top, hh, sub=sub)
    top = next_top(top + max(hb, hm))
    _, h = stepping(m, PAD, top, SPAN)
    return top + h + PAD


def main() -> None:
    global drawing
    manifest = json.loads((PANELS / "panels.json").read_text(encoding="utf-8"))
    drawing = Drawing("f25", PANELS, paths.VAULT)
    height = layout(manifest)
    drawing.check()
    drawing.write(OUT)
    logger.info("canvas %.0f x %.0f units -> %s", WIDTH, height, OUT.name)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
