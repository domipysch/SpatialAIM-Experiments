"""
Panels of the Results 2.4 figure -- biological and topological plausibility --
each as its own PNG, for assembly in Excalidraw
(00_Obsidian/Excalidraw/figure_2_4.excalidraw.md).

Same division of labour as figures 2.1-2.3: the images carry only what is
anchored to data. Letters, headings and legends are drawn in Excalidraw, and
everything else belongs in the caption.

Two rows of three, one row per criterion; every cell carries its own letter:

    a  sketch_spatial     what the criterion means, drawn on a toy tissue patch
    b  plaus_nhood_...    SpatialAIM against the annotation, neighbourhood enrichment
    c  plaus_local_...    the same, local spatial purity
    d  sketch_coherence   the same idea on a toy expression graph
    e  plaus_modularity_st_expression       modularity over the slice's spots
    f  plaus_modularity_shared              modularity over the reference's cells

    dk_bar               the colour bar all four scatters share

The two sketches use one device: the *same* points and the *same* label counts,
once as observed and once permuted. That is exactly the null the numbers are
measured against -- the two spatial metrics are z-scores over that shuffle, and
modularity is the excess over a degree-preserving random partition -- so the
sketches define the criterion and show its chance level in the same breath.
Their only difference is where the nodes sit: a tissue rectangle in row (a), an
expression embedding in row (b). Same question, two graphs.

The scatters use **one mapper**, `MAPPER` below. Both axes of a panel always run
the same aligner -- SpatialAIM's proposal against the annotation's identity cut -- so
the mapper is a paired factor and never a confound; showing all five turned each
panel into 50 points for 10 datasets and read as clutter. That the result holds
for every aligner belongs in the text, where the numbers are.

The two spatial panels are logarithmic on both axes: their z-scores span two
decades across references, and on a log axis the offset from the diagonal is a
ratio. The two modularity panels are linear -- they live inside a single decade,
where a log axis buys nothing. Chance is 0 on all four and therefore off-scale
below every panel; every run is above it, which the caption states.

Aggregation is per reference, never pooled over pairs: 71 of the 92 pairs are
retina islands.

Usage
-----
    python sections/04_biology/panels.py                    # the figure
    python sections/04_biology/panels.py --mapper tangram   # the same figure,
                                                            # that mapper
    python sections/04_biology/panels.py --all-mappers      # a trial

Another ``--mapper`` writes a complete copy of the figure into its own folder,
``figure_2_4_<mapper>``, and never touches the panels the figure itself embeds.
``--all-mappers`` stays a trial: all five at once is 50 points per panel, so it
goes to the preview folder under suffixed names and writes no manifest.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import explore as B                           # noqa: E402  (sibling script)
from spatialaim_experiments import exp1 as L
from spatialaim_experiments import paths
from spatialaim_experiments.style import (
    CMAP_DK,
    DK_LIMIT,
    INK,
    MAPPERS,
    PALE,
    PanelWriter,
    use_paper_rcparams,
    bare,
)

logger = logging.getLogger(__name__)

FIGURE = "2_4"
PREVIEW_ROOT = paths.outputs("04_biology")

# The mapper the figure shows. One aligner, ten points per panel -- see the
# module docstring.
MAPPER = "nearest_centroid"

# Figures 2.3 and 2.4 are authored at the same point sizes -- this figure used
# to inherit them by importing the reconstruction panel module, which applied
# them as a side effect of being imported.
use_paper_rcparams()

# set per invocation by ``main``; these are the figure's own
panels = PanelWriter(paths.panel_dir(FIGURE), PREVIEW_ROOT)


def folders(mapper: str) -> tuple[Path, Path]:
    """Where a mapper's panels go. The figure's mapper keeps the unsuffixed
    folders; every other one gets its own copy beside them, so a copy can never
    overwrite a panel the figure embeds."""
    if mapper == MAPPER:
        return paths.panel_dir(FIGURE), PREVIEW_ROOT
    return paths.panel_dir(FIGURE, mapper), PREVIEW_ROOT / mapper

# The scatter axes is square -- both axes carry the same measure on the same
# scale, so anything else would bend the parity diagonal off 45 degrees. The
# rect is derived from the margins rather than typed, so ``set_aspect("equal")``
# has nothing left to shrink and the four panels tile exactly.
SIDE = 1.58                       # print inches of the plotting box itself
# The right margin holds the second half of the last x tick label: the axis now
# ends *on* a tick, so "1000" overhangs the box by half its width. Left and
# bottom each hold a two-line axis title -- the measure over the side of the
# comparison -- which is what the box gave up 0.07 in of its side for. The top
# holds nothing: every y title is now shorter, set rotated, than the box it is
# centred on. It was 0.18 while panel b's read "neighbourhood enrichment  z",
# which overran the box and had its "z" cut off by the edge of the image.
MARGIN = dict(left=0.66, bottom=0.53, right=0.18, top=0.08)
SCATTER_SIZE = (SIDE + MARGIN["left"] + MARGIN["right"],
                SIDE + MARGIN["bottom"] + MARGIN["top"])
SCATTER_RECT = (MARGIN["left"] / SCATTER_SIZE[0], MARGIN["bottom"] / SCATTER_SIZE[1],
                SIDE / SCATTER_SIZE[0], SIDE / SCATTER_SIZE[1])
# Wide enough for the label under it, not for the gradient. 3.45 rather than the
# 2.95 it was: "SpatialAIM" is five characters longer than the name it replaced,
# and the label was being clipped at both ends.
BAR_SIZE = (3.45, 0.66)

# The sketch column is *shorter* than a scatter, not equal: it needs a foot for
# its two corner arrows, where the scatter needs one for tick labels plus an axis
# title. Both start at the row's top edge, so their plotting bands still coincide
# -- the sketch simply stops earlier.
SKETCH_SIZE = (2.05, 0.0)   # filled in below, once the axis margins are known

# The three cell-type hues of figures 2.1 and 2.2, reused in the sketches so a
# reader who has seen those recognises what the dots stand for.
SKETCH_COLOURS = ("#D55E00", "#0072B2", "#009E73")

# Print inches reserved at the left and foot of a sketch for its corner arrows.
# The left needs more: "UMAP 2" is set rotated, so its line height is what has to
# clear the vertical arrow.
AXIS_LEFT, AXIS_FOOT = 0.30, 0.26
FS_AXIS = 8.0
SKETCH_HALF_LABEL = 0.16          # the strip "observed" / "chance" is set on
SKETCH_SIZE = (SKETCH_SIZE[0], MARGIN["top"] + SIDE + AXIS_FOOT)
# The arrows say which space the nodes are laid out in -- tissue in row (a), an
# embedding in row (b) -- which is the only difference between the two sketches
# and the thing a reader would otherwise have to infer. They are a schematic
# frame, not a measured one: the layouts are synthetic, and modularity is a
# property of the graph's edges, not of any 2-D projection of them. The caption
# has to say so.
SKETCH_AXES = {"sketch_spatial": ("x", "y"),
               "sketch_coherence": ("UMAP 1", "UMAP 2")}

# Cell-type difference: CMAP_DK and DK_LIMIT come from spatialaim_experiments.style,
# the same scale and clip figure 2.3 uses for the same quantity -- purple = SpatialAIM
# used fewer cell types, green = more.

# The two lines an axis is titled with, per measure, as (x, y). The panel's
# letter names the measure too, in the drawing, but an axis title is what makes a
# scatter readable on its own -- and it has to be short enough to sit against a
# 1.58 in box. The two z-scores say "z" so nobody reads them as the raw quantity;
# the two modularities name the graph they are measured on, since that is exactly
# what separates them.


def _titles(measure: str) -> tuple[str, str]:
    """The measure over the side of the comparison, on both axes."""
    return f"{measure}\nwhen mapping annotation", f"{measure}\nwhen using SpatialAIM"


AXIS = {
    # Panel b is the user's own wording. Its caption spells the measure out
    # directly above the panel, so the axes carry the short form, and the side of
    # the comparison is what is read first going up the y axis.
    "nhood_mean_self_zscore": _titles("mean z-score"),
    "local_purity_zscore": _titles("z-score"),
    "modularity_st_expression": _titles("modularity"),
    "modularity_shared": _titles("modularity"),
}

def save(fig, name: str, preview_only: bool = False) -> None:
    """Write a panel, with one escape hatch.

    ``preview_only`` points the write at the preview folder instead, so a trial
    rendering cannot land in the Excalidraw folder and silently replace a panel
    the drawing embeds."""
    writer = PanelWriter(panels.preview) if preview_only else panels
    writer.save(fig, name)


# ---------------------------------------------------------------------------
# (a) SpatialAIM against the annotation
# ---------------------------------------------------------------------------

def cells(df: pd.DataFrame, metric: str, mapper: str | None = None) -> pd.DataFrame:
    """One row per (reference, mapper): the reference's mean over its slices, for
    SpatialAIM, for the annotation baseline, and for the difference in cell types.

    ``mapper`` keeps only that one, leaving ten points per panel."""
    pr = B.frames(df, metric)
    if mapper is not None:
        pr = pr[pr["mapper"] == mapper]
    return (pr.groupby(["reference", "mapper"])
              .agg(spatialaim=(metric, "mean"), ann=(f"{metric}_base", "mean"),
                   dk=("dk", "mean"))
              .reset_index())


# Only the two spatial z-scores span decades and need a logarithmic axis. The
# two modularities live inside a single decade, where a log axis buys nothing and
# costs the reader a linear reading of a bounded quantity.
LOG_METRICS = {"nhood_mean_self_zscore", "local_purity_zscore"}


# The limits sit this fraction of the axis beyond the outermost tick. Without it
# a point whose value is close to a snapped tick has its marker cut in half by
# the spine -- the top point of the neighbourhood-enrichment panel sat 0.02
# decades below its 300 tick and lost its upper third.
AXIS_PAD = 0.03


def _linear_axis(lo: float, hi: float) -> tuple[tuple[float, float], list[float]]:
    """Limits snapped outward to a round step, with a tick on every step.

    The ticks are equidistant by construction and the data's extremes sit inside
    the first and last of them, so the axis carries no tick-free tail and no
    label lands in the corner. The limits then clear the outermost ticks by
    ``AXIS_PAD`` so a marker sitting on one is drawn whole."""
    ticks = mpl.ticker.MaxNLocator(nbins=6, steps=[1, 2, 2.5, 5, 10]
                                   ).tick_values(lo, hi)
    values = [float(v) for v in ticks]
    pad = (values[-1] - values[0]) * AXIS_PAD
    return (values[0] - pad, values[-1] + pad), values


def _log_axis(lo: float, hi: float) -> tuple[tuple[float, float], list[float]]:
    """The same, on a logarithmic axis: limits snapped outward to half-decades,
    with a tick on every one.

    Half-decades are labelled 1, 3, 10, 30 rather than 1, 3.16, 10, 31.6, so the
    rungs sit 0.477 and 0.523 decades apart instead of a uniform 0.5. That is 5%
    of one gap -- invisible on the page -- and it buys labels a reader can read.
    Snapping outward also replaces the old multiplicative padding, so the ends of
    the axis are ticks rather than arbitrary values."""
    ladder = [m * 10.0 ** e
              for e in range(int(np.floor(np.log10(lo))) - 1,
                             int(np.ceil(np.log10(hi))) + 1)
              for m in (1.0, 3.0)]
    start = max(v for v in ladder if v <= lo)
    end = min(v for v in ladder if v >= hi)
    pad = (end / start) ** AXIS_PAD
    return (start / pad, end * pad), [v for v in ladder if start <= v <= end]


def panel_scatter(cell: pd.DataFrame, metric: str, name: str,
                  preview_only: bool = False) -> None:
    """Both axes are titled with the measure over the side of the comparison, so
    a panel says what it plots without its letter. The letter and its caption are
    drawn in Excalidraw, above the panel, which keeps all four images the same
    geometry so they tile without adjustment."""
    fig = plt.figure(figsize=SCATTER_SIZE)
    ax = fig.add_axes(list(SCATTER_RECT))

    lo = float(min(cell["spatialaim"].min(), cell["ann"].min()))
    hi = float(max(cell["spatialaim"].max(), cell["ann"].max()))
    log = metric in LOG_METRICS
    span, ticks = (_log_axis if log else _linear_axis)(lo, hi)
    ax.plot(span, span, color=PALE, lw=0.9, zorder=1)      # parity

    ax.scatter(cell["ann"], cell["spatialaim"], c=np.clip(cell["dk"], -DK_LIMIT, DK_LIMIT),
               cmap=CMAP_DK, vmin=-DK_LIMIT, vmax=DK_LIMIT, s=34,
               linewidths=0.5, edgecolors=INK, zorder=2)

    if log:                            # set the scale first: it resets locators
        ax.set_xscale("log")
        ax.set_yscale("log")
    # x and y carry the same quantity on the same scale, so they get the same
    # limits and the same ticks; the parity diagonal then runs corner to corner.
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(mpl.ticker.FixedLocator(ticks))
        axis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        axis.set_minor_locator(mpl.ticker.NullLocator())
    ax.set_xlim(*span)
    ax.set_ylim(*span)
    ax.set_aspect("equal")
    x_title, y_title = AXIS[metric]
    ax.set_xlabel(x_title, labelpad=1.5, linespacing=1.05)
    ax.set_ylabel(y_title, labelpad=2.0, linespacing=1.05)
    ax.tick_params(pad=1.5)
    bare(ax)
    save(fig, name, preview_only)


# ---------------------------------------------------------------------------
# the two sketches: what each criterion means, and what chance looks like
# ---------------------------------------------------------------------------

def _knn_edges(points: np.ndarray, k: int) -> np.ndarray:
    """Each point's k nearest neighbours, as unique undirected pairs."""
    d = np.linalg.norm(points[:, None, :] - points[None, :, :], axis=-1)
    np.fill_diagonal(d, np.inf)
    neighbours = np.argsort(d, axis=1)[:, :k]
    pairs = {(min(i, j), max(i, j))
             for i, row in enumerate(neighbours) for j in row}
    return np.array(sorted(pairs))


def _draw_graph(ax, points: np.ndarray, labels: np.ndarray,
                edges: np.ndarray) -> None:
    """Nodes coloured by label, edges within a label in that label's colour and
    edges across labels pale -- so a coherent labelling shows as coloured
    filaments and an incoherent one as a pale mesh."""
    for i, j in edges:
        same = labels[i] == labels[j]
        ax.plot(points[[i, j], 0], points[[i, j], 1],
                color=SKETCH_COLOURS[labels[i]] if same else "#E4E4E4",
                lw=0.75 if same else 0.55, zorder=1 if same else 0,
                solid_capstyle="round")
    ax.scatter(points[:, 0], points[:, 1],
               c=[SKETCH_COLOURS[v] for v in labels], s=11,
               linewidths=0.35, edgecolors="white", zorder=2)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.margins(0.03)


def _corner_arrows(ax, x_label: str, y_label: str) -> None:
    """Two short arrows outside the lower-left corner of ``ax``, naming the space
    the nodes are laid out in. Drawn once per sketch, since both halves are the
    same points."""
    # Offsets are given in print inches and converted, because the two sketches'
    # boxes are not the same shape and an offset in axes fractions drifts with
    # them -- which is how the label ended up drawn on top of the arrow.
    box = ax.get_position()
    fig_w, fig_h = ax.figure.get_size_inches()
    w_in, h_in = box.width * fig_w, box.height * fig_h

    arrow = dict(arrowstyle="-|>,head_width=0.13,head_length=0.30",
                 lw=0.7, color=INK, shrinkA=0, shrinkB=0)
    common = dict(xycoords="axes fraction", textcoords="axes fraction",
                  annotation_clip=False)
    ax.annotate("", xy=(0.26, -0.06 / h_in), xytext=(0.0, -0.06 / h_in),
                arrowprops=arrow, **common)
    ax.annotate("", xy=(-0.06 / w_in, 0.62), xytext=(-0.06 / w_in, 0.0),
                arrowprops=arrow, **common)
    ax.text(0.13, -0.13 / h_in, x_label, transform=ax.transAxes,
            fontsize=FS_AXIS, ha="center", va="top", color=INK)
    ax.text(-0.15 / w_in, 0.31, y_label, transform=ax.transAxes,
            fontsize=FS_AXIS, ha="center", va="center", rotation=90,
            color=INK)


def _sketch_rects(label_in: float = SKETCH_HALF_LABEL,
                  ) -> list[tuple[float, float, float, float]]:
    """The two halves' axes rects, stacked inside the vertical band a scatter's
    plotting box occupies, each with room for its name above it.

    Matching that band is what makes a row read as one thing: the sketch and the
    two scatters then start and end at the same height on the page. The lower
    half gives up ``AXIS_FOOT`` at the bottom and both give up ``AXIS_LEFT`` at
    the left, for the two corner arrows -- the halves are the same points in the same
    coordinates, so one arrow pair at the foot of the column serves both."""
    height = SKETCH_SIZE[1]
    left = AXIS_LEFT / SKETCH_SIZE[0]
    label = label_in / height
    foot = AXIS_FOOT / height
    patch = (SIDE / height - 2 * label) / 2       # SIDE = a scatter's box, in
    width = 0.98 - left                           # inches, so the bands coincide
    return [(left, foot + patch + label, width, patch),
            (left, foot, width, patch)]


def _sketch(points: np.ndarray, labels: np.ndarray, edges: np.ndarray,
            name: str, preview_only: bool = False) -> list[float]:
    """One sketch: the same points and the same label counts twice, observed
    above and permuted below.

    Returns the top edge of each half as a fraction of the image height from the
    top. The two halves are named in the drawing, on the strip of empty image
    each rect leaves above it -- a label to the left would need more width than
    the narrowest column of the figure has."""
    fig = plt.figure(figsize=SKETCH_SIZE)
    rects = _sketch_rects()
    shuffled = np.random.default_rng(7).permutation(labels)
    axes = [fig.add_axes(list(rect)) for rect in rects]
    for ax, values in zip(axes, (labels, shuffled)):
        _draw_graph(ax, points, values, edges)
    _corner_arrows(axes[-1], *SKETCH_AXES[name])
    save(fig, name, preview_only)
    return [1 - (bottom + height) for _x, bottom, _w, height in rects]


def sketch_spatial(preview_only: bool = False) -> list[float]:
    """Spatial organisation, on a toy tissue patch: a jittered grid of cells in
    three contiguous territories, joined by the spatial neighbour graph. Under
    the permutation the territories dissolve into salt and pepper -- the same
    picture Results 2.2 shows for the invisible split."""
    rng = np.random.default_rng(0)
    nx, ny = 19, 7            # 2.7:1, the aspect the stacked rects leave
    grid = np.stack(np.meshgrid(np.arange(nx), np.arange(ny)), -1).reshape(-1, 2)
    points = grid + rng.uniform(-0.32, 0.32, grid.shape)
    # three territories as the Voronoi cells of three seeds, so they are
    # contiguous by construction rather than by luck
    seeds = np.array([[3.0, 4.6], [9.5, 1.7], [15.5, 4.6]])
    labels = np.argmin(np.linalg.norm(points[:, None] - seeds[None], axis=-1), axis=1)
    return _sketch(points, labels, _knn_edges(points, 4), "sketch_spatial",
                   preview_only)


def sketch_coherence(preview_only: bool = False) -> list[float]:
    """Transcriptional coherence, on a toy expression embedding: three clusters
    of cells joined by their expression k-NN graph. Same device, same colours --
    only the layout changes, which is the whole point of the row pair."""
    rng = np.random.default_rng(3)
    centres = np.array([[0.0, 0.0], [5.4, 1.1], [10.4, -0.5]])
    points = np.concatenate([
        centre + rng.normal(0, 1.0, (34, 2)) for centre in centres])
    labels = np.repeat(np.arange(3), 34)
    return _sketch(points, labels, _knn_edges(points, 3), "sketch_coherence",
                   preview_only)


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
    for spine in bar.spines.values():
        spine.set_visible(False)
    bar.tick_params(length=1.5, pad=1.0)
    save(fig, "dk_bar")


# ---------------------------------------------------------------------------

def main(mapper: str | None = MAPPER) -> None:
    """A named ``mapper`` writes that mapper's figure, into its own folder unless
    it is the figure's own. ``None`` keeps all five, which is a trial: it goes to
    the preview folder under suffixed names and writes no manifest."""
    global panels
    df = L.load()
    trial = mapper is None
    suffix = "__all" if trial else ""
    panels = PanelWriter(*folders(mapper or MAPPER))
    stats = {}
    for metric in B.METRICS:
        cell = cells(df, metric, mapper)
        panel_scatter(cell, metric, f"plaus_{metric}{suffix}", preview_only=trial)
        pr = B.frames(df, metric)
        if mapper is not None:
            pr = pr[pr["mapper"] == mapper]
        stats[metric] = {
            "above_parity": int((cell["spatialaim"] > cell["ann"]).sum()),
            "cells": int(len(cell)),
            "above_chance_runs": int((pr[metric] > 0).sum()),
            "runs": int(len(pr)),
            "rho_delta_dk": float(pr["delta"].corr(pr["dk"], method="spearman")),
            "dk_min": float(cell["dk"].min()),
            "dk_max": float(cell["dk"].max()),
        }
    for metric, s in stats.items():
        logger.info("%-26s above parity %2d/%d   rho(delta, dk) %+.2f",
                    metric, s["above_parity"], s["cells"], s["rho_delta_dk"])
    if trial:
        logger.info("trial for all mappers -> %s (the figure is untouched)",
                    panels.preview)
        return

    tops = {"sketch_spatial": sketch_spatial(),
            "sketch_coherence": sketch_coherence()}
    panel_dk_bar()
    panels.write_manifest({
        "mapper": mapper,
        "limits": {"dk": DK_LIMIT},
        "stats": stats,
        # top edge of each sketch half, as a fraction of the image height from
        # the top, so the drawing can set "observed"/"chance" just above it
        "sketch_tops": tops,
        "sketch_labels": ["observed", "chance"],
    })


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    ap = argparse.ArgumentParser(description=__doc__)
    group = ap.add_mutually_exclusive_group()
    group.add_argument("--mapper", choices=MAPPERS, default=MAPPER,
                       help=f"which mapper's ten points to draw (default {MAPPER}, "
                            "the figure itself); any other value writes a copy of "
                            "the figure into figure_2_4_<mapper>")
    group.add_argument("--all-mappers", action="store_true",
                       help="trial: all five mappers, 50 points per panel")
    args = ap.parse_args()
    main(None if args.all_mappers else args.mapper)
