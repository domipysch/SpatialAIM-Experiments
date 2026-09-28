"""
Panels of the biology/topology figure (figure 4) for one mapper, each as its own
PNG, for assembly in Excalidraw (Excalidraw/figure_4/<mapper>/figure_4_<mapper>.excalidraw.md).

Two rows of three, one row per criterion:

    a  sketch_spatial                   what spatial organisation means, on a toy tissue patch
    b  plaus_nhood_mean_self_zscore     SpatialAIM against the annotation, neighbourhood enrichment
    c  plaus_local_purity_zscore        the same, local spatial purity
    d  sketch_coherence                 what transcriptional coherence means, on a toy expression graph
    e  plaus_modularity_st_expression   modularity over the slice's spots
    f  plaus_modularity_shared          modularity over the reference's cells
       dk_bar                           the colour bar all four scatters share
"""

from __future__ import annotations
import argparse
import logging
import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from spatialaim_experiments import exp1, paths
from spatialaim_experiments.style import (
    CMAP_DK,
    DK_LIMIT,
    INK,
    MAPPERS,
    PALE,
    PanelWriter,
    bare,
    use_paper_rcparams,
)

logger = logging.getLogger(__name__)

FIGURE = "4"
TREE = "average"

# The four metrics, each with the measure its axis titles name. Only the two
# spatial z-scores span decades and get logarithmic axes; the modularities live
# inside one decade.
METRICS = {
    "nhood_mean_self_zscore": "mean z-score",
    "local_purity_zscore": "z-score",
    "modularity_st_expression": "modularity",
    "modularity_shared": "modularity",
}
LOG_METRICS = {"nhood_mean_self_zscore", "local_purity_zscore"}

use_paper_rcparams()
# set by main() to the mapper's panel folder, figure_4/<mapper>/panels/
panels: PanelWriter

# Scatter: a square plotting box (x and y are the same quantity, so parity runs at
# 45 degrees), with margins for the two-line axis titles (left, bottom) and the
# overhang of the last x tick label (right).
SIDE = 1.58
MARGIN = dict(left=0.66, bottom=0.53, right=0.18, top=0.08)
SCATTER_SIZE = (SIDE + MARGIN["left"] + MARGIN["right"], SIDE + MARGIN["bottom"] + MARGIN["top"])
SCATTER_AXES = (MARGIN["left"] / SCATTER_SIZE[0], MARGIN["bottom"] / SCATTER_SIZE[1],
                SIDE / SCATTER_SIZE[0], SIDE / SCATTER_SIZE[1])
# The limits sit this fraction of the axis beyond the outermost tick, so a point
# on a tick is not cut in half by the spine.
AXIS_PAD = 0.03

# Sketch: the same vertical band as a scatter's plotting box (so a row reads as
# one), plus a foot and a left margin for the two corner arrows, and a strip above
# each half for its name ("observed" / "chance", set in Excalidraw).
AXIS_LEFT, AXIS_FOOT = 0.30, 0.26
HALF_LABEL = 0.16
SKETCH_SIZE = (2.05, MARGIN["top"] + SIDE + AXIS_FOOT)
FS_AXIS = 8.0
# The three cell-type hues of figures 1 and 2, so the dots read as cell types.
SKETCH_COLOURS = ("#D55E00", "#0072B2", "#009E73")
# The space the nodes are laid out in -- schematic, the layouts are synthetic.
SKETCH_AXES = {"sketch_spatial": ("x", "y"), "sketch_coherence": ("UMAP 1", "UMAP 2")}

BAR_SIZE = (3.45, 0.66)


# ---------------------------------------------------------------------------
# the four scatters: SpatialAIM against the annotation
# ---------------------------------------------------------------------------

def comparisons(df: pd.DataFrame, metric: str, mapper: str) -> pd.DataFrame:
    """Each SpatialAIM run of ``mapper`` next to its annotation baseline.

    Input
        df      best_per_run.csv (exp1.load()): one row per run
        metric  e.g. "local_purity_zscore"
        mapper  e.g. "wann"

    Output
        one row per pair (average tree, ``mapper`` only), with the columns of
        exp1.paired() -- ``metric`` and ``k`` for SpatialAIM, ``<metric>_base``
        and ``k_base`` for the annotation, ``delta`` -- plus
          dk   k - k_base, how many more cell types SpatialAIM used
    """
    out = exp1.paired(df, metric)
    out = out[(out["tree"] == TREE) & (out["mapper"] == mapper)].reset_index(drop=True)
    out["dk"] = out["k"] - out["k_base"]
    return out


def scatter_points(pr: pd.DataFrame, metric: str) -> pd.DataFrame:
    """The points of one scatter panel.

    Input
        pr      the output of comparisons(): one row per pair
        metric  the same metric comparisons() was called with

    Output
        one row per reference (one point in the panel), each averaged over the
        reference's pairs:
          reference    the reference's name
          spatialaim   ``metric`` for SpatialAIM          (the point's y)
          ann          ``metric`` for the annotation      (the point's x)
          dk           the difference in cell types       (the point's colour)
    """
    return (pr.groupby("reference")
              .agg(spatialaim=(metric, "mean"), ann=(f"{metric}_base", "mean"), dk=("dk", "mean"))
              .reset_index())


def linear_axis(lo: float, hi: float) -> tuple[tuple[float, float], list[float]]:
    """Limits snapped outward to a round step, a tick on every step."""
    ticks = [float(v) for v in
             mpl.ticker.MaxNLocator(nbins=6, steps=[1, 2, 2.5, 5, 10]).tick_values(lo, hi)]
    pad = (ticks[-1] - ticks[0]) * AXIS_PAD
    return (ticks[0] - pad, ticks[-1] + pad), ticks


def log_axis(lo: float, hi: float) -> tuple[tuple[float, float], list[float]]:
    """Limits snapped outward to half-decades labelled 1, 3, 10, 30, ..., a tick
    on every one."""
    ladder = [m * 10.0 ** e
              for e in range(int(np.floor(np.log10(lo))) - 1, int(np.ceil(np.log10(hi))) + 1)
              for m in (1.0, 3.0)]
    start = max(v for v in ladder if v <= lo)
    end = min(v for v in ladder if v >= hi)
    pad = (end / start) ** AXIS_PAD
    return (start / pad, end * pad), [v for v in ladder if start <= v <= end]


def panel_scatter(cell: pd.DataFrame, metric: str) -> None:
    """x = annotation, y = SpatialAIM, same limits and ticks on both, colour = dk."""
    fig = plt.figure(figsize=SCATTER_SIZE)
    ax = fig.add_axes(SCATTER_AXES)
    lo = float(min(cell["spatialaim"].min(), cell["ann"].min()))
    hi = float(max(cell["spatialaim"].max(), cell["ann"].max()))
    log = metric in LOG_METRICS
    span, ticks = (log_axis if log else linear_axis)(lo, hi)
    ax.plot(span, span, color=PALE, lw=0.9, zorder=1)      # parity
    ax.scatter(cell["ann"], cell["spatialaim"], c=np.clip(cell["dk"], -DK_LIMIT, DK_LIMIT),
               cmap=CMAP_DK, vmin=-DK_LIMIT, vmax=DK_LIMIT, s=34,
               linewidths=0.5, edgecolors=INK, zorder=2)
    if log:                            # set the scale first: it resets the locators
        ax.set_xscale("log")
        ax.set_yscale("log")
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(mpl.ticker.FixedLocator(ticks))
        axis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        axis.set_minor_locator(mpl.ticker.NullLocator())
    ax.set_xlim(*span)
    ax.set_ylim(*span)
    ax.set_aspect("equal")
    measure = METRICS[metric]
    ax.set_xlabel(f"{measure}\nwhen mapping annotation", labelpad=1.5, linespacing=1.05)
    ax.set_ylabel(f"{measure}\nwhen using SpatialAIM", labelpad=2.0, linespacing=1.05)
    ax.tick_params(pad=1.5)
    bare(ax)
    panels.save(fig, f"plaus_{metric}")


def panel_dk_bar() -> None:
    """The one colour bar the four scatters share."""
    fig = plt.figure(figsize=BAR_SIZE)
    bar = fig.add_axes([0.16, 0.58, 0.68, 0.24])
    bar.imshow(np.linspace(-DK_LIMIT, DK_LIMIT, 256)[None, :], cmap=CMAP_DK,
               aspect="auto", extent=(-DK_LIMIT, DK_LIMIT, 0, 1))
    bar.set_yticks([])
    bar.set_xticks([-DK_LIMIT, 0, DK_LIMIT])
    bar.set_xticklabels([f"≤{-DK_LIMIT}", "0", f"≥+{DK_LIMIT}"])
    bar.set_xlabel("Δ number of cell types  (SpatialAIM − annotation)", labelpad=2.5)
    bar.spines[:].set_visible(False)
    bar.tick_params(length=1.5, pad=1.0)
    panels.save(fig, "dk_bar")


# ---------------------------------------------------------------------------
# the two sketches: what each criterion means, and what chance looks like
# ---------------------------------------------------------------------------

def knn_edges(points: np.ndarray, k: int) -> np.ndarray:
    """Each point's k nearest neighbours, as unique undirected pairs."""
    d = np.linalg.norm(points[:, None, :] - points[None, :, :], axis=-1)
    np.fill_diagonal(d, np.inf)
    neighbours = np.argsort(d, axis=1)[:, :k]
    return np.array(sorted({(min(i, j), max(i, j))
                            for i, row in enumerate(neighbours) for j in row}))


def draw_graph(ax, points: np.ndarray, labels: np.ndarray, edges: np.ndarray) -> None:
    """Nodes coloured by label; edges within a label in its colour, edges across
    labels pale -- a coherent labelling shows as coloured filaments."""
    for i, j in edges:
        same = labels[i] == labels[j]
        ax.plot(points[[i, j], 0], points[[i, j], 1],
                color=SKETCH_COLOURS[labels[i]] if same else "#E4E4E4",
                lw=0.75 if same else 0.55, zorder=1 if same else 0, solid_capstyle="round")
    ax.scatter(points[:, 0], points[:, 1], c=[SKETCH_COLOURS[v] for v in labels], s=11,
               linewidths=0.35, edgecolors="white", zorder=2)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[:].set_visible(False)
    ax.margins(0.03)


def corner_arrows(ax, x_label: str, y_label: str) -> None:
    """Two short arrows outside the lower-left corner of ``ax``, naming the space
    the nodes are laid out in. Offsets are in print inches, converted to axes
    fractions, since the two sketches' boxes differ in shape."""
    box = ax.get_position()
    fig_w, fig_h = ax.figure.get_size_inches()
    w_in, h_in = box.width * fig_w, box.height * fig_h
    arrow = dict(arrowstyle="-|>,head_width=0.13,head_length=0.30",
                 lw=0.7, color=INK, shrinkA=0, shrinkB=0)
    common = dict(xycoords="axes fraction", textcoords="axes fraction", annotation_clip=False)
    ax.annotate("", xy=(0.26, -0.06 / h_in), xytext=(0.0, -0.06 / h_in), arrowprops=arrow, **common)
    ax.annotate("", xy=(-0.06 / w_in, 0.62), xytext=(-0.06 / w_in, 0.0), arrowprops=arrow, **common)
    ax.text(0.13, -0.13 / h_in, x_label, transform=ax.transAxes,
            fontsize=FS_AXIS, ha="center", va="top", color=INK)
    ax.text(-0.15 / w_in, 0.31, y_label, transform=ax.transAxes,
            fontsize=FS_AXIS, ha="center", va="center", rotation=90, color=INK)


def sketch_rects() -> list[tuple[float, float, float, float]]:
    """The axes of the two halves (observed above, permuted below), stacked in the
    band of a scatter's plotting box, each with its name strip above it."""
    height = SKETCH_SIZE[1]
    left = AXIS_LEFT / SKETCH_SIZE[0]
    label = HALF_LABEL / height
    foot = AXIS_FOOT / height
    half = (SIDE / height - 2 * label) / 2
    width = 0.98 - left
    return [(left, foot + half + label, width, half), (left, foot, width, half)]


def sketch(points: np.ndarray, labels: np.ndarray, edges: np.ndarray, name: str) -> list[float]:
    """One sketch: the same graph with the observed labels above and permuted
    labels below. Returns the top edge of each half as a fraction of the image
    height from the top, where make_excalidraw.py sets the half's name."""
    fig = plt.figure(figsize=SKETCH_SIZE)
    rects = sketch_rects()
    axes = [fig.add_axes(rect) for rect in rects]
    shuffled = np.random.default_rng(7).permutation(labels)
    for ax, values in zip(axes, (labels, shuffled)):
        draw_graph(ax, points, values, edges)
    corner_arrows(axes[-1], *SKETCH_AXES[name])
    panels.save(fig, name)
    return [1 - (bottom + height) for _, bottom, _, height in rects]


def sketch_spatial() -> list[float]:
    """A toy tissue patch: a jittered grid in three contiguous territories (the
    Voronoi cells of three seeds), joined by its spatial neighbour graph."""
    rng = np.random.default_rng(0)
    grid = np.stack(np.meshgrid(np.arange(19), np.arange(7)), -1).reshape(-1, 2)
    points = grid + rng.uniform(-0.32, 0.32, grid.shape)
    seeds = np.array([[3.0, 4.6], [9.5, 1.7], [15.5, 4.6]])
    labels = np.argmin(np.linalg.norm(points[:, None] - seeds[None], axis=-1), axis=1)
    return sketch(points, labels, knn_edges(points, 4), "sketch_spatial")


def sketch_coherence() -> list[float]:
    """A toy expression embedding: three clusters of cells joined by their k-NN graph."""
    rng = np.random.default_rng(3)
    centres = np.array([[0.0, 0.0], [5.4, 1.1], [10.4, -0.5]])
    points = np.concatenate([centre + rng.normal(0, 1.0, (34, 2)) for centre in centres])
    labels = np.repeat(np.arange(3), 34)
    return sketch(points, labels, knn_edges(points, 3), "sketch_coherence")


# ---------------------------------------------------------------------------

def main(mapper: str) -> None:
    """All panels of ``mapper``'s figure, into figure_4/<mapper>/panels/."""
    global panels
    panels = PanelWriter(paths.panel_dir(FIGURE, mapper))
    df = exp1.load()
    for metric in METRICS:
        pr = comparisons(df, metric, mapper)
        cell = scatter_points(pr, metric)
        panel_scatter(cell, metric)
        logger.info("%-26s SpatialAIM above the annotation for %d/%d references; "
                    "above chance in %d/%d runs", metric,
                    (cell["spatialaim"] > cell["ann"]).sum(), len(cell), (pr[metric] > 0).sum(), len(pr))
    panel_dk_bar()
    panels.write_manifest({"sketch_tops": {"sketch_spatial": sketch_spatial(),
                                           "sketch_coherence": sketch_coherence()}})


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mapper", required=True, choices=MAPPERS,
                    help="whose figure to draw, into figure_4/<mapper>/")
    main(ap.parse_args().mapper)
