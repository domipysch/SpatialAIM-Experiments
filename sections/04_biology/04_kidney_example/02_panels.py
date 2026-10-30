"""
Step 2 of figure 4's kidney example: what spatial assortativity looks like, as
PNGs for the drawing.

Reads the two mappings 01_run_mappings.py computed (annotation, SpatialAIM) and
draws, per side:

    map_<side>    the whole section, small, every spot in the colour of its state;
                  the zoom window is boxed
    zoom_<side>   the window, large, with the neighbourhood graph drawn: an edge
                  joining two spots of the same state in that state's colour, an
                  edge between two states pale grey. Assortativity is exactly the
                  share of these same-state edges, corrected for chance.

The graph is SpatialAIM's own (squidpy spatial kNN, K_SPATIAL neighbours), and the
window's assortativity (global and mean per state) is computed with SpatialAIM's
own function on the graph rebuilt from the window's spots.

State colours follow the annotated types: every type its own colour; a SpatialAIM state takes the colour of the type dominating its reference
cells, a second state dominated by the same type a turned hue.

The images go to the vault, next to the nearest-centroid drawing of figure 4:

    Excalidraw/figure_4/<mapper>/kidney_example/   *.png + panels.json

``panels.json`` also carries what the drawing writes around them: each side's K and
its metrics on the whole section and on the window. No state is named: the
figure is about spatial assortativity, not about which cell types these are.

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
from matplotlib.collections import LineCollection
from matplotlib.colors import to_hex, to_rgb

from spatialaim.io import start_folder_name
from spatialaim.metrics.topology import assortativity

from spatialaim_experiments import paths
from spatialaim_experiments.style import MAPPERS, PanelWriter, use_paper_rcparams

logger = logging.getLogger(__name__)
use_paper_rcparams()

FIGURE = "4"
RESULTS = paths.RESULTS / "biology_example"

# Print width of one map, in inches; the height follows the section's aspect.
# The maps are small -- they only locate the zoom; the zoom carries the figure.
MAP_WIDTH = 1.0
# Marker area in pt^2: 126k spots on 1.0 in.
SPOT_SIZE = 0.02

# The zoom pair: one window of the section, the same in both maps, at spot level,
# where a mixed and a coherent mapping differ visibly whatever the palette.
# (x0, x1, y0, y1) in the ST coordinates, per PairID; a pair without a window
# gets no zoom. The kidney's spans medulla, the NAPSA-high proximal-tubule ring and cortex.
ZOOM_WINDOWS = {4: (41000.0, 50000.0, 38000.0, 45875.0)}
ZOOM_WIDTH = 1.95     # inches; ~4k spots
# The spatial graph: SpatialAIM's analysis builds a kNN graph with this many
# neighbours (spatialaim.analysis.topology.K_SPATIAL).
K_SPATIAL = 6
# Zoom edges: same-state edges in the state's colour, the others pale.
EDGE_SAME = 0.4         # line width, pt
EDGE_CROSS = 0.3
EDGE_PALE = "#DCDCDC"
ZOOM_BOX = dict(edgecolor="#1e1e1e", linewidth=0.6)

# Colours for the annotated types, most abundant type first. All saturated and
# of similar weight: a pale colour (sky blue, yellow) makes a map look washed out
# next to one without it, which is the palette talking, not the mapping.
# 22 entries, one per type of the largest kidney annotation (the 22 TMS classes),
# so no two types share a colour; the last eleven are the less distinct ones and go
# to the rarest types.
TYPE_COLOURS = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00",
                "#B2182B", "#7B3294", "#8C6D31", "#4D4D4D", "#999933", "#1B7837",
                "#E7298A", "#66A61E", "#00A6A6", "#A6761D", "#3F48CC", "#F0027F",
                "#8E0152", "#2D6A4F", "#B35806", "#5E3C99", "#C51B7D"]
# The second (third, ...) state of the same dominant type turns the type's hue by
# this much (of a full turn) at the same luminance -- a darker shade instead
# would hand one map a near-black the other does not have.
HUE_TURNS = [0.0, 0.08, -0.08]

# Shorter names for the kidney's labels, so the key fits beside a map; a type not
# listed keeps its own name. The 22 TMS classes (cell_ontology_class), then the
# targets only Liu et al.'s relabeling (cell_ontology_relabel) adds.
SHORT_NAMES = {
    "B cell": "B cell",
    "NK cell": "NK cell",
    "T cell": "T cell",
    "brush cell": "Brush cell",
    "epithelial cell of proximal tubule": "Proximal tubule",
    "fenestrated cell": "Fenestrated",
    "fibroblast": "Fibroblast",
    "kidney capillary endothelial cell": "Capillary endothelial",
    "kidney cell": "Kidney cell",
    "kidney collecting duct epithelial cell": "Collecting duct",
    "kidney collecting duct principal cell": "Collecting duct principal",
    "kidney cortex artery cell": "Cortex artery",
    "kidney distal convoluted tubule epithelial cell": "Distal convoluted tubule",
    "kidney loop of Henle ascending limb epithelial cell": "Henle ascending limb",
    "kidney loop of Henle thick ascending limb epithelial cell": "Henle thick asc. limb",
    "kidney mesangial cell": "Mesangial",
    "kidney proximal convoluted tubule epithelial cell": "Proximal convoluted",
    "kidney proximal straight tubule epithelial cell": "Proximal straight",
    "leukocyte": "Leukocyte",
    "lymphocyte": "Lymphocyte",
    "macrophage": "Macrophage",
    "podocyte": "Podocyte",
    "glomerular capillary endothelial cell": "Glomerular capillary",
    "kidney blood vessel cell": "Blood vessel",
    "kidney loop of Henle epithelial cell": "Loop of Henle",
    "pericyte": "Pericyte",
    "plasma cell": "Plasma cell",
}
# A type contributes to a state's name in the key if it makes up at least this
# share of the state's reference cells.
NAME_SHARE = 0.15


# ---------------------------------------------------------------------------
# reading one mapping
# ---------------------------------------------------------------------------

def run_root(folder: Path, side: str, ann_key: str | None, mapper: str) -> Path:
    """One side's run root as 01_run_mappings.py wrote it:
    ``<side>/<start>/<mapper>``, start ``leiden`` or ``annotation_<key>``."""
    return folder / side / ("leiden" if ann_key is None else start_folder_name(ann_key)) / mapper


def load_mapping(root: Path, k: int) -> tuple[pd.Series, np.ndarray]:
    """One run's result at K.

    Output
        cell_state   per reference cell (index = cell id) the state it belongs to
        spot_state   per ST spot the state it is mapped to (argmax of P)
    """
    # the start clustering is shared by the start's mappers, one folder up
    start = ad.read_h5ad(root.parent / "start_clustering.h5ad", backed="r").obs["start_cluster"]
    cut = pd.read_csv(root / f"k_{k:03d}" / "start_cluster_to_state.csv")
    cell_state = start.astype(int).map(dict(zip(cut["start_cluster"], cut["state"])))
    P = ad.read_h5ad(root / f"k_{k:03d}" / "spot_to_state_mapping_soft.h5ad")
    return cell_state, np.asarray(P.X).argmax(axis=1)


def metrics(root: Path, k: int) -> dict:
    """The figure's four metrics for this mapping, from its k_comparison.csv."""
    row = pd.read_csv(root / "k_comparison.csv").set_index("k").loc[k]
    return {m: float(row[m]) for m in ["spatial_assort_celltype", "spatial_assort_global",
                                       "modularity_st_expression", "modularity_shared"]}


def knn_graph(xy: np.ndarray):
    """SpatialAIM's spatial graph over these spots: squidpy kNN, K_SPATIAL
    neighbours (its connectivities, as the analysis step uses them)."""
    import squidpy as sq  # heavy; only here

    a = ad.AnnData(np.zeros((len(xy), 1), dtype=np.float32))
    a.obsm["spatial"] = xy
    sq.gr.spatial_neighbors_knn(a, n_neighs=K_SPATIAL)
    return a.obsp["spatial_connectivities"].tocoo()


def window_metrics(labels: np.ndarray, graph) -> dict:
    """Global and mean per-state assortativity on the window, with SpatialAIM's
    own function (graph built from the window's spots), named as in
    k_comparison.csv."""
    a = assortativity(labels, graph)
    return {"spatial_assort_global": a.value, "spatial_assort_celltype": a.mean_state}


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

def section_axes(xy: np.ndarray, width: float, window: tuple | None = None):
    """A frameless figure of ``width`` inches over the spots' extent (or over
    ``window``), y downwards as in the GUI."""
    if window is None:
        (x0, y0), (x1, y1) = xy.min(axis=0), xy.max(axis=0)
    else:
        x0, x1, y0, y1 = window
    fig = plt.figure(figsize=(width, width * (y1 - y0) / (x1 - x0)))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)
    ax.set_aspect("equal")
    return fig, ax


def map_panel(panels: PanelWriter, name: str, xy: np.ndarray, colours: list[str],
              window: tuple | None) -> None:
    """The section at MAP_WIDTH inches, every spot in its state's colour, drawn in
    one random order (fixed seed) so no state lies on top by its id; the zoom
    ``window`` outlined."""
    fig, ax = section_axes(xy, MAP_WIDTH)
    order = np.random.default_rng(0).permutation(len(xy))
    ax.scatter(xy[order, 0], xy[order, 1], s=SPOT_SIZE, linewidths=0, rasterized=True,
               c=[colours[i] for i in order])
    if window is not None:
        wx0, wx1, wy0, wy1 = window
        ax.add_patch(plt.Rectangle((wx0, wy0), wx1 - wx0, wy1 - wy0, fill=False,
                                   zorder=3, **ZOOM_BOX))
    ax.set_axis_off()
    panels.save(fig, name)


def zoom_panel(panels: PanelWriter, name: str, xy: np.ndarray, labels: np.ndarray,
               colour_of: dict, graph, window: tuple) -> None:
    """The window at ZOOM_WIDTH inches with its neighbourhood graph: the edges
    between two states pale and underneath, the same-state edges on top in their
    state's colour; framed like its box on the map."""
    fig, ax = section_axes(xy, ZOOM_WIDTH, window)
    segments = np.stack([xy[graph.row], xy[graph.col]], axis=1)
    same = labels[graph.row] == labels[graph.col]
    ax.add_collection(LineCollection(segments[~same], colors=EDGE_PALE,
                                     linewidths=EDGE_CROSS, rasterized=True))
    ax.add_collection(LineCollection(segments[same], linewidths=EDGE_SAME, rasterized=True,
                                     colors=[colour_of[v] for v in labels[graph.row][same]]))
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_edgecolor(ZOOM_BOX["edgecolor"])
        spine.set_linewidth(ZOOM_BOX["linewidth"])
    panels.save(fig, name)


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
    ann_root = run_root(folder, "annotation", ann["ann_key"], args.mapper)
    ann_cells, ann_spots = load_mapping(ann_root, ann["k"])
    ann_states = describe_states(ann_cells, cell_type, colours)
    ann_colour = {s: colours[t] for s, t in ann_states["dominant"].items()}

    # SpatialAIM
    aim_root = run_root(folder, "spatialaim", aim["ann_key"], args.mapper)
    aim_cells, aim_spots = load_mapping(aim_root, aim["k"])
    aim_states = describe_states(aim_cells, cell_type, colours)
    logger.info("SpatialAIM states:\n%s", aim_states.to_string())

    panels = PanelWriter(paths.figure_dir(FIGURE, args.mapper) / "kidney_example")
    window = ZOOM_WINDOWS.get(args.pair)
    if window is not None:
        wx0, wx1, wy0, wy1 = window
        inside = ((xy[:, 0] >= wx0) & (xy[:, 0] <= wx1)
                  & (xy[:, 1] >= wy0) & (xy[:, 1] <= wy1))
        window_graph = knn_graph(xy[inside])

    manifest = {"selection": selection, "zoom": window}
    for side, spec, root, labels, colour_of in (
            ("annotation", ann, ann_root, ann_spots, ann_colour),
            ("spatialaim", aim, aim_root, aim_spots, aim_states["colour"].to_dict())):
        map_panel(panels, f"map_{side}", xy, [colour_of[v] for v in labels], window)
        manifest[side] = {"k": spec["k"], "metrics": metrics(root, spec["k"])}
        if window is not None:
            zoom_panel(panels, f"zoom_{side}", xy[inside], labels[inside], colour_of,
                       window_graph, window)
            manifest[side]["window_metrics"] = window_metrics(labels[inside], window_graph)
        logger.info("%-10s section assortativity %.3f | window %s", side,
                    manifest[side]["metrics"]["spatial_assort_global"],
                    manifest[side].get("window_metrics"))
    panels.write_manifest(manifest)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
