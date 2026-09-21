"""
Panels of the Results 2.3 reconstruction figure, each as its own PNG, for
assembly in Excalidraw (00_Obsidian/Excalidraw/figure_2_3.excalidraw.md).

Same approach as the robustness figure: the panels carry only what is anchored to data
(cells, axes, colour bars). Letters, headings, the mapper names and the legends
are drawn in Excalidraw; everything else belongs in the caption.

    a  delta_spot      AIM (proposed K) minus annotation (K = T), per reference
       delta_gene      x mapper x tree; the columns are references
       delta_k         and the difference in how many cell types that took
    b  mappers         the five mappers at their own proposal, centred within
                       each reference so the references can be compared at all

Row (a) is one argument in three panels: reconstruction barely moves, and AIM
gets there with substantially fewer cell types.

Absolute cosine similarity is not comparable across references (33-gene panels
against 981-gene panels, and provenance sets a per-pair ceiling), so every panel
is a paired comparison *within* a reference. Aggregation is per reference, never
pooled over pairs: 71 of the 92 pairs are retina islands.

Sizes are in final print inches at the figure's width in the paper; panels.json
records what was written so make_excalidraw.py can place each image at its true
size.

Usage
-----
    python sections/03_reconstruction/panels.py
"""

from __future__ import annotations

import logging

import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from aim_experiments import exp1 as L
from aim_experiments import paths
from aim_experiments.style import (
    CMAP_DK,
    DK_LIMIT,
    MAPPER_COLOURS,
    MAPPER_LABELS,
    PALE,
    REFERENCE_LABELS as SHORT,
    TREE_LABELS,
    PanelWriter,
    use_paper_rcparams,
    bare,
)

logger = logging.getLogger(__name__)

FIGURE = "2_3"
OUT = paths.panel_dir(FIGURE)
PREVIEW = paths.outputs("03_reconstruction")

# Row (b) groups by configuration (metric x linkage) and colours by mapper, so
# the five mappers are read against each other inside one group and the same
# five colours recur in every group -- see aim_experiments.style for the palette
# itself, which every figure of the paper shares.

# The figure prints at the paper's 5.17 in text width, so the drawing is scaled
# by 0.71 and every point size is authored 1/0.71 larger than it should print.
use_paper_rcparams()

METRIC_SPOT = "cossim_hard_raw_spot"
METRIC_GENE = "cossim_hard_raw_gene"

# Three heat maps of identical geometry in one row. The mapper names are drawn
# once, in Excalidraw, to the left of the row -- putting them inside the first
# image would make its cells narrower than the other two's, and then the three
# panels could not be read against each other.
# Wider than it was: row (a)'s left column no longer spells the mapper out, only
# a coloured dot and a linkage label, so the width that used to hold "nearest
# centroid" goes back to the heat maps. All three grow equally, so their cells
# stay aligned.
HEAT_SIZE = (2.13, 2.65)
# A rotated y-label is as tall as its longest line is long: "vs mean across
# mappers" needs 1.55 in, so the plotting box has to be at least that or the
# label is silently clipped at both ends.
MAPPER_SIZE = (7.0, 2.14)
MAPPER_RECT = (0.105, 0.215, 0.885, 0.757)

# Reconstruction deltas are shown on a fixed scale rather than one fitted to the
# data: the point of the panel is that almost every cell is close to zero, and
# an autoscaled bar would magnify that near-zero noise into strong colour and
# say the opposite. 0.25 sits just above the largest cell in either panel
# (per spot -0.134 .. +0.240, per gene -0.075 .. +0.078), so nothing is clipped
# and the outliers still stand out against a mostly pale field.
DELTA_LIMIT = 0.25
# DK_LIMIT (the cell-type difference, clipped at 20) is shared with figure 2.4
# and lives in aim_experiments.style.

panels = PanelWriter(OUT, PREVIEW)


# ---------------------------------------------------------------------------
# (a) AIM against the annotation, over the whole design space
# ---------------------------------------------------------------------------

def _matrix(df: pd.DataFrame, metric: str, value: str) -> pd.DataFrame:
    """One cell per (mapper x tree, reference): the reference's mean over its
    slices. ``value`` is 'delta' (the difference in ``metric``) or 'dk' (the
    difference in how many cell types it took).

    Both are AIM minus the annotation, in that order -- the same subtraction in
    every panel of the row, so one sign convention covers all three and a
    negative cell always means the annotation side was larger."""
    pr = L.paired(df, metric)
    pr["dk"] = pr["k"] - pr["k_base"]
    pr["row"] = pr["mapper"] + " · " + pr["tree"]
    cell = pr.groupby(["row", "reference"])[value].mean().unstack("reference")
    rows = [f"{m} · {t}" for m in L.MAPPERS for t in L.TREES]
    return cell.reindex(index=[r for r in rows if r in cell.index],
                        columns=[c for c in SHORT if c in cell.columns])


def _delta_matrix(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    return _matrix(df, metric, "delta")


# blue = AIM ahead, orange = the annotation ahead; the figure's own two families.
#
# The scale is deliberately **not** linear in colour position. Its stops are
# anchored in cosine-similarity units, and the near-white band is stretched to
# +-0.1: almost every cell of row (a) lies inside that, which is the row's whole
# point, and a linear ramp turned that agreement into visible colour and said the
# opposite. Past 0.1 the ramp saturates quickly, so the handful of references
# that really do differ are the only ones carrying colour. The colour bar is
# drawn from the same map over a linear axis and therefore shows the stretch
# itself, with a tick at +-0.1 marking where the pale band ends.
# Crimson and teal, not the orange and blue this scale used to have: the mappers
# are named in row (a) by a coloured dot now, and those dots are the figure's
# oranges, blues, gold and green. A heat map in the same two hues sat right
# beside them and the two readings kept being confused.
#
# The figure has five categorical mapper hues plus the purple/green cell-type
# scale, so no hue pair is entirely free. Crimson is the closest thing to it, and
# it reads as "the annotation is ahead" more directly than a magenta does. Its
# one weakness: at full saturation crimson sits near Tangram's vermillion, so the
# few strongly negative cells and the Tangram dot are similar reds.
DELTA_STOPS = [
    #  |delta|   annotation ahead   AIM ahead
    (0.000, "#FFFFFF", "#FFFFFF"),
    (0.050, "#FCF0F1", "#EBF7F5"),
    (0.100, "#F4C7CE", "#C7E6E1"),
    (0.175, "#B2182B", "#35978F"),
    (0.250, "#67001F", "#01665E"),
]


def _diverging(name: str, stops: list[tuple[float, str, str]], limit: float):
    """A diverging colormap from stops given in data units, symmetric about 0.

    ``stops`` ascends from 0; each entry is (|value|, colour below zero, colour
    above zero). A value at the limit lands on the end of the map, so changing
    ``limit`` rescales the whole scale rather than silently clipping the anchors.
    """
    nodes = [(0.5 - value / (2 * limit), negative)
             for value, negative, _positive in reversed(stops)]
    nodes += [(0.5 + value / (2 * limit), positive)
              for value, _negative, positive in stops[1:]]
    return mpl.colors.LinearSegmentedColormap.from_list(name, nodes)


CMAP_DELTA = _diverging("aim_delta", DELTA_STOPS, DELTA_LIMIT)

# The third scale, for the cell-type difference, is CMAP_DK in
# aim_experiments.style: a third hue pair, so nobody reads cell types as a
# cosine, and shared with figure 2.4, which plots the same quantity.

# Shared by all three heat maps so their cells align. The top quarter is room
# for the rotated reference names -- "colorectal" is the longest of them and was
# being clipped -- and the bottom fifth holds the colour bar, its ticks and the
# two-line label naming what the colours actually are.
AXES_RECT = (0.02, 0.265, 0.96, 0.475)
BAR_RECT = (0.12, 0.185, 0.76, 0.028)


def panel_heat(matrix: pd.DataFrame, name: str, limit: float, cmap,
               ticks: list, tick_labels: list[str], unit: str) -> list[float]:
    """One of the three heat maps of row (a). Geometry is identical across the
    three; the mapper names and linkage labels live in the drawing, not in the
    image.

    Returns where every row is centred, as a fraction of the image height from
    the top, so make_excalidraw.py can hang a label off each one and off each
    pair's midpoint.
    """
    fig = plt.figure(figsize=HEAT_SIZE)
    left, bottom, width, height = AXES_RECT
    ax = fig.add_axes([left, bottom, width, height])
    ax.imshow(np.clip(matrix.to_numpy(), -limit, limit), cmap=cmap,
              vmin=-limit, vmax=limit, aspect="auto", interpolation="nearest")

    # the reference names go above the matrix: below, they would collide with
    # the colour bar, and rotated labels need room a caption line cannot give
    ax.xaxis.set_ticks_position("top")
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_xticklabels([SHORT[c] for c in matrix.columns], rotation=90,
                       ha="center", va="bottom")
    ax.set_yticks([])
    ax.tick_params(length=0, pad=1.5)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for y in range(2, matrix.shape[0], 2):     # separate the five mappers
        ax.axhline(y - 0.5, color="white", lw=1.6)

    # inset, so the labels on the two end ticks stay inside the image
    bar = fig.add_axes(list(BAR_RECT))
    bar.imshow(np.linspace(-limit, limit, 256)[None, :], cmap=cmap, aspect="auto",
               extent=(-limit, limit, 0, 1))
    bar.set_yticks([])
    bar.set_xticks(ticks)
    bar.set_xticklabels(tick_labels)
    # The bar carries its own unit, so the row above needs only one heading and
    # the reader never has to match a panel to a caption clause to know what the
    # colours mean.
    bar.set_xlabel(unit, labelpad=2.5)
    for spine in bar.spines.values():
        spine.set_visible(False)
    bar.tick_params(length=1.5, pad=1.0)
    panels.save(fig, name)

    n = matrix.shape[0]
    return [1 - bottom - height * (n - i - 0.5) / n for i in range(n)]


# ---------------------------------------------------------------------------
# (d) the five mappers at their own proposal
# ---------------------------------------------------------------------------

def panel_mappers(df: pd.DataFrame) -> None:
    """Within each reference the five mappers are centred on their common mean,
    so what is plotted is how a mapper does *relative to the others on the same
    data* -- the only way five mappers can be compared across ten references
    whose absolute cosine similarities are not commensurable."""
    aim = df[df["role"] == "aim"].copy()

    def centred(metric: str) -> pd.DataFrame:
        """One value per (reference, mapper, tree): how far that mapper sits from
        the mean of all five on the same pair and linkage, averaged over the
        reference's slices. Ten references -> ten values per box."""
        per_pair = (aim.groupby(["reference", "pair_dir", "tree", "mapper"])[metric]
                       .mean().reset_index())
        per_pair["centred"] = per_pair[metric] - per_pair.groupby(
            ["reference", "pair_dir", "tree"])[metric].transform("mean")
        return (per_pair.groupby(["reference", "mapper", "tree"])["centred"]
                        .mean().reset_index())

    # Four groups, one per configuration; five boxes in each, one per mapper.
    groups = [(centred(metric), metric_label, tree)
              for metric, metric_label in ((METRIC_SPOT, "per spot"),
                                           (METRIC_GENE, "per gene"))
              for tree in L.TREES]

    fig = plt.figure(figsize=MAPPER_SIZE)
    ax = fig.add_axes(list(MAPPER_RECT))
    ax.axhline(0, color=PALE, lw=0.8, zorder=1)

    slot = 0.15
    step = 0.17          # five boxes span 0.68, leaving 0.32 between groups
    for g, (cell, _label, tree) in enumerate(groups):
        for m, mapper in enumerate(L.MAPPERS):
            values = cell[(cell["mapper"] == mapper)
                          & (cell["tree"] == tree)]["centred"].to_numpy()
            colour = MAPPER_COLOURS[mapper]
            box = ax.boxplot([values], positions=[g + (m - 2) * step], widths=slot,
                             patch_artist=True, showfliers=False,
                             medianprops=dict(color=colour, lw=1.2),
                             whiskerprops=dict(color=colour, lw=0.7),
                             capprops=dict(color=colour, lw=0.7))
            patch = box["boxes"][0]
            patch.set_facecolor(colour)
            patch.set_alpha(0.5)
            patch.set_edgecolor(colour)
            patch.set_linewidth(0.9)

    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([f"{label}\n{TREE_LABELS[tree]} linkage"
                        for _cell, label, tree in groups])
    ax.set_xlim(-0.55, len(groups) - 0.45)
    ax.set_ylabel("Δ cosine similarity\nvs mean across mappers")
    bare(ax)
    panels.save(fig, "mappers")


# ---------------------------------------------------------------------------

def main() -> None:
    df = L.load()
    # The middle pair of ticks marks where the pale band ends, so the stretched
    # scale can be read off the bar instead of taken on trust. No "<=" on the
    # ends: unlike the cell-type panel, this scale clips nothing.
    pale = DELTA_STOPS[2][0]
    delta_ticks = [-DELTA_LIMIT, -pale, 0, pale, DELTA_LIMIT]
    delta_labels = [f"{-DELTA_LIMIT:g}", f"{-pale:g}", "0",
                    f"+{pale:g}", f"+{DELTA_LIMIT:g}"]
    rows_y = panel_heat(_delta_matrix(df, METRIC_SPOT), "delta_spot",
                        DELTA_LIMIT, CMAP_DELTA, delta_ticks, delta_labels,
                        "Δ cosine similarity per spot\n(AIM − annotation)")
    panel_heat(_delta_matrix(df, METRIC_GENE), "delta_gene",
               DELTA_LIMIT, CMAP_DELTA, delta_ticks, delta_labels,
               "Δ cosine similarity per gene\n(AIM − annotation)")
    dk = _matrix(df, METRIC_SPOT, "dk")
    panel_heat(dk, "delta_k", DK_LIMIT, CMAP_DK,
               [-DK_LIMIT, 0, DK_LIMIT],
               [f"≤{-DK_LIMIT}", "0", f"≥+{DK_LIMIT}"],
               "Δ number of cell types\n(AIM − annotation)")
    panel_mappers(df)

    panels.write_manifest({
        "limits": {"delta": DELTA_LIMIT, "dk": DK_LIMIT,
                   "dk_max": float(np.nanmax(dk.to_numpy())),
                   "dk_min": float(np.nanmin(dk.to_numpy()))},
        # where each heat-map row is centred, as a fraction of the image height
        # from the top: every row for the linkage labels, and each mapper pair's
        # midpoint for the mapper name
        "heat_row_lines": rows_y,
        "heat_rows": [(a + b) / 2 for a, b in zip(rows_y[::2], rows_y[1::2])],
        "trees": [TREE_LABELS[t] for t in L.TREES],
        "legends": {
            "mappers": [{"key": m, "label": MAPPER_LABELS[m],
                         "colour": MAPPER_COLOURS[m]} for m in L.MAPPERS],
        },
    })


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
