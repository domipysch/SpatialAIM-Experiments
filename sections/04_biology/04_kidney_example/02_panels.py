"""
Step 2 of figure 4's kidney example: the two spatial maps, as PNGs for the drawing.

Reads the two mappings 01_run_mappings.py computed and draws the ST section twice,
coloured by the annotation (left) and by SpatialAIM's states (right):

    map_annotation   every spot coloured by the annotated type it is mapped to
    map_spatialaim   every spot coloured by its SpatialAIM state
    zoom_annotation  one window of the section (ZOOM_WINDOWS), outlined on both maps,
    zoom_spatialaim  at spot level

Both maps share one set of colours, so they can be compared directly:

  * every annotated type has its own colour, assigned by how many reference cells
    the type has (most abundant first)
  * a SpatialAIM state takes the colour of the annotated type that dominates its
    reference cells -- a state merging several types shows as one colour where
    the annotation shows several; a second state dominated by the same type gets
    a darker shade of it

Spots are drawn in one random order (fixed seed), not state by state, so where
spots overlap the colour on top reflects how many spots each state has there.

The images go to the vault, next to the nearest-centroid drawing of figure 4:

    Excalidraw/figure_4/<mapper>/kidney_example/   map_*.png + panels.json

``panels.json`` also carries what the drawing writes around the maps: each side's
K, its key (label and colour per type or state, with its share of spots) and its
metrics.

    python sections/04_biology/04_kidney_example/02_panels.py [--pair 4]
"""

from __future__ import annotations

import argparse
import colorsys
import json
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import to_hex, to_rgb

from spatialaim_experiments import paths
from spatialaim_experiments.style import MAPPERS, PanelWriter, use_paper_rcparams

logger = logging.getLogger(__name__)
use_paper_rcparams()

FIGURE = "4"
RESULTS = paths.RESULTS / "biology_example"

# Print width of one map, in inches; the height follows the section's aspect.
MAP_WIDTH = 2.3
# Marker area in pt^2. 126k spots on 2.3 in: large enough to fill the tissue,
# small enough that neighbouring spots of different states stay apart.
SPOT_SIZE = 0.06

# The zoom pair: one window of the section, the same in both maps, at spot level,
# where a mixed and a coherent mapping differ visibly whatever the palette.
# (x0, x1, y0, y1) in the ST coordinates, per PairID; a pair without a window
# gets no zoom. The kidney's crosses the cortex-medulla border.
ZOOM_WINDOWS = {4: (14000.0, 30000.0, 26000.0, 40000.0)}
ZOOM_WIDTH = 1.15     # inches, half a map
ZOOM_SPOT = 1.0       # marker area in pt^2: ~7k spots on 1.15 in
ZOOM_BOX = dict(edgecolor="#1e1e1e", linewidth=0.6)

# Colours for the annotated types, most abundant type first. All saturated and
# of similar weight: a pale colour (sky blue, yellow) makes a map look washed out
# next to one without it, which is the palette talking, not the mapping.
TYPE_COLOURS = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00",
                "#B2182B", "#7B3294", "#8C6D31", "#4D4D4D", "#999933", "#1B7837"]
# The second (third, ...) state of the same dominant type turns the type's hue by
# this much (of a full turn) at the same luminance -- a darker shade instead
# would hand one map a near-black the other does not have.
HUE_TURNS = [0.0, 0.08, -0.08]

# Shorter names for the kidney's cell-ontology labels, so the key fits beside a
# map; a type not listed keeps its own name.
SHORT_NAMES = {
    "epithelial cell of proximal tubule": "Proximal tubule",
    "kidney loop of Henle epithelial cell": "Loop of Henle",
    "endothelial cell": "Endothelial",
    "immune cell": "Immune",
    "kidney collecting duct epithelial cell": "Collecting duct",
    "kidney distal convoluted tubule epithelial cell": "Distal convoluted tubule",
    "podocyte": "Podocyte",
    "pericyte": "Pericyte",
    "fibroblast": "Fibroblast",
}
# A type contributes to a state's name in the key if it makes up at least this
# share of the state's reference cells.
NAME_SHARE = 0.15


# ---------------------------------------------------------------------------
# reading one mapping
# ---------------------------------------------------------------------------

def load_mapping(root: Path, k: int) -> tuple[pd.Series, np.ndarray]:
    """One run's result at K.

    Output
        cell_state   per reference cell (index = cell id) the state it belongs to
        spot_state   per ST spot the state it is mapped to (argmax of P)
    """
    start = ad.read_h5ad(root / "start_clustering.h5ad", backed="r").obs["start_cluster"]
    cut = pd.read_csv(root / f"k_{k:03d}" / "start_cluster_to_state.csv")
    cell_state = start.astype(int).map(dict(zip(cut["start_cluster"], cut["state"])))
    P = ad.read_h5ad(root / f"k_{k:03d}" / "spot_to_state_mapping_soft.h5ad")
    return cell_state, np.asarray(P.X).argmax(axis=1)


def metrics(root: Path, k: int) -> dict:
    """The figure's four metrics for this mapping, from its k_comparison.csv."""
    row = pd.read_csv(root / "k_comparison.csv").set_index("k").loc[k]
    return {m: float(row[m]) for m in ["nhood_mean_self_zscore", "local_purity_zscore",
                                       "modularity_st_expression", "modularity_shared"]}


# ---------------------------------------------------------------------------
# one set of colours for both maps
# ---------------------------------------------------------------------------

def luminance(colour: str) -> float:
    """Relative luminance (WCAG): how light a colour looks, 0 = black, 1 = white."""
    def channel(v: float) -> float:
        return v / 12.92 if v <= 0.03928 else ((v + 0.055) / 1.055) ** 2.4
    r, g, b = map(channel, to_rgb(colour))
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def variant(colour: str, turn: float) -> str:
    """``colour`` with its hue turned by ``turn`` (of a full turn), its lightness
    then set so the luminance matches the original's (0 = unchanged)."""
    if turn == 0:
        return colour
    h, _l, s = colorsys.rgb_to_hls(*to_rgb(colour))
    h = (h + turn) % 1.0
    target, lo, hi = luminance(colour), 0.0, 1.0
    for _ in range(40):                       # luminance rises with lightness
        mid = (lo + hi) / 2
        if luminance(to_hex(colorsys.hls_to_rgb(h, mid, s))) < target:
            lo = mid
        else:
            hi = mid
    return to_hex(colorsys.hls_to_rgb(h, (lo + hi) / 2, s))


def type_colours(cell_type: pd.Series) -> dict[str, str]:
    """A colour per annotated type, most abundant type first."""
    order = cell_type.value_counts().index
    return {t: TYPE_COLOURS[i % len(TYPE_COLOURS)] for i, t in enumerate(order)}


def describe_states(cell_state: pd.Series, cell_type: pd.Series,
                    colours: dict[str, str]) -> pd.DataFrame:
    """Each state's dominant annotated type, colour and name for the key.

    Output
        one row per state: dominant, colour, name (its types of at least
        NAME_SHARE, largest first, joined by " + "), cells
    """
    table = pd.crosstab(cell_state, cell_type.loc[cell_state.index])
    share = table.div(table.sum(axis=1), axis=0)
    out = pd.DataFrame({"dominant": share.idxmax(axis=1), "cells": table.sum(axis=1)})
    out["name"] = [" + ".join(t for t, v in share.loc[s].sort_values(ascending=False).items()
                              if v >= NAME_SHARE)
                   for s in out.index]
    # states sharing a dominant type: the largest keeps its colour, the others
    # turn its hue
    out["colour"] = ""
    for dominant, group in out.groupby("dominant"):
        for turn, state in zip(HUE_TURNS, group.sort_values("cells", ascending=False).index):
            out.loc[state, "colour"] = variant(colours[dominant], turn)
    # states sharing a name are numbered, largest first, so the key tells them apart
    for _name, group in out.groupby("name"):
        if len(group) > 1:
            for i, state in enumerate(group.sort_values("cells", ascending=False).index, 1):
                out.loc[state, "name"] = f"{out.loc[state, 'name']} ({i})"
    return out


# ---------------------------------------------------------------------------
# drawing
# ---------------------------------------------------------------------------

def spots(ax, xy: np.ndarray, colours: list[str], size: float) -> None:
    """All spots in one random order (fixed seed), so no state is drawn on top
    of the others by its id."""
    order = np.random.default_rng(0).permutation(len(xy))
    ax.scatter(xy[order, 0], xy[order, 1], s=size, linewidths=0, rasterized=True,
               c=[colours[i] for i in order])


def map_panel(panels: PanelWriter, name: str, xy: np.ndarray, colours: list[str],
              window: tuple | None) -> None:
    """The section, frameless, at MAP_WIDTH inches; y downwards as in the GUI.
    With a zoom ``window`` its outline is drawn on the section."""
    (x0, y0), (x1, y1) = xy.min(axis=0), xy.max(axis=0)
    fig = plt.figure(figsize=(MAP_WIDTH, MAP_WIDTH * (y1 - y0) / (x1 - x0)))
    ax = fig.add_axes([0, 0, 1, 1])
    spots(ax, xy, colours, SPOT_SIZE)
    if window is not None:
        wx0, wx1, wy0, wy1 = window
        ax.add_patch(plt.Rectangle((wx0, wy0), wx1 - wx0, wy1 - wy0, fill=False,
                                   zorder=3, **ZOOM_BOX))
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)
    ax.set_aspect("equal")
    ax.set_axis_off()
    panels.save(fig, name)


def zoom_panel(panels: PanelWriter, name: str, xy: np.ndarray, colours: list[str],
               window: tuple) -> None:
    """The zoom ``window`` at ZOOM_WIDTH inches, framed like its box on the map."""
    wx0, wx1, wy0, wy1 = window
    inside = (xy[:, 0] >= wx0) & (xy[:, 0] <= wx1) & (xy[:, 1] >= wy0) & (xy[:, 1] <= wy1)
    fig = plt.figure(figsize=(ZOOM_WIDTH, ZOOM_WIDTH * (wy1 - wy0) / (wx1 - wx0)))
    ax = fig.add_axes([0, 0, 1, 1])
    spots(ax, xy[inside], [c for c, keep in zip(colours, inside) if keep], ZOOM_SPOT)
    ax.set_xlim(wx0, wx1)
    ax.set_ylim(wy1, wy0)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_edgecolor(ZOOM_BOX["edgecolor"])
        spine.set_linewidth(ZOOM_BOX["linewidth"])
    panels.save(fig, name)


def key(labels: np.ndarray, colour_of: dict, name_of: dict) -> list[dict]:
    """The key of one map: per type or state its label, colour and share of
    spots, largest share first."""
    share = pd.Series(labels).value_counts(normalize=True)
    return [{"label": name_of[v], "colour": colour_of[v], "share": round(float(s), 4)}
            for v, s in share.items()]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pair", type=int, default=4, help="PairID (default: 4, the kidney)")
    ap.add_argument("--mapper", choices=MAPPERS, default="nearest_centroid")
    ap.add_argument("--tree", choices=["average", "ward"], default="average")
    args = ap.parse_args()

    folders = list(RESULTS.glob(f"{args.pair:03d}_*/{args.mapper}_{args.tree}"))
    if not folders:
        raise SystemExit("no mappings found - run 01_run_mappings.py first")
    folder = folders[0]
    selection = json.loads((folder / "selection.json").read_text(encoding="utf-8"))
    ann, aim = selection["annotation"], selection["spatialaim"]

    st = ad.read_h5ad(paths.DATASETS / "ST" / f"{selection['st']}.h5ad", backed="r")
    xy = np.asarray(st.obsm["spatial"], dtype=float)
    sc = ad.read_h5ad(paths.DATASETS / "scRNA" / f"{selection['reference']}.h5ad", backed="r")
    cell_type = sc.obs[ann["ann_key"]].astype(str).map(lambda t: SHORT_NAMES.get(t, t))
    colours = type_colours(cell_type)

    # the annotation: state = annotated type (identity cut at K = T)
    ann_root = folder / "annotation" / args.mapper
    ann_cells, ann_spots = load_mapping(ann_root, ann["k"])
    ann_states = describe_states(ann_cells, cell_type, colours)
    ann_colour = {s: colours[t] for s, t in ann_states["dominant"].items()}
    ann_name = ann_states["dominant"].to_dict()

    # SpatialAIM
    aim_root = folder / "spatialaim" / args.mapper
    aim_cells, aim_spots = load_mapping(aim_root, aim["k"])
    aim_states = describe_states(aim_cells, cell_type, colours)
    logger.info("SpatialAIM states:\n%s", aim_states.to_string())

    panels = PanelWriter(paths.figure_dir(FIGURE, args.mapper) / "kidney_example")
    window = ZOOM_WINDOWS.get(args.pair)
    for side, colours_ in (("annotation", [ann_colour[s] for s in ann_spots]),
                           ("spatialaim", [aim_states.loc[s, "colour"] for s in aim_spots])):
        map_panel(panels, f"map_{side}", xy, colours_, window)
        if window is not None:
            zoom_panel(panels, f"zoom_{side}", xy, colours_, window)

    panels.write_manifest({
        "selection": selection,
        "zoom": window,
        "annotation": {"k": ann["k"], "metrics": metrics(ann_root, ann["k"]),
                       "key": key(ann_spots, ann_colour, ann_name)},
        "spatialaim": {"k": aim["k"], "metrics": metrics(aim_root, aim["k"]),
                       "key": key(aim_spots, aim_states["colour"].to_dict(),
                                  aim_states["name"].to_dict())},
    })


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
