"""
Panels of the reconstruction figure (figure 3), each as its own PNG, for assembly
in Excalidraw (Excalidraw/figure_3/figure_3.excalidraw.md).

    a  delta_spot   Δ cosine similarity per spot, SpatialAIM (proposed K) - annotation (K = T)
       delta_gene   the same per gene
       delta_k      Δ number of cell types it took
                    rows: mapper x tree; columns: the ten references
    b  mappers      the five mappers at their own proposal, centred within each pair

Absolute cosine similarity is not comparable across references (33-gene against
981-gene panels), so every panel is a comparison *within* a reference, averaged
per reference - never pooled over pairs.
"""

from __future__ import annotations
import logging
import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from spatialaim_experiments import exp1, paths
from spatialaim_experiments.style import (
    CMAP_DK,
    DK_LIMIT,
    HEAT_AXES,
    MAPPER_COLOURS,
    MAPPERS,
    PALE,
    REFERENCE_LABELS,
    TREE_LABELS,
    TREES,
    PanelWriter,
    bare,
    use_paper_rcparams,
)

FIGURE = "3"
METRIC_SPOT, METRIC_GENE = "cossim_hard_raw_spot", "cossim_hard_raw_gene"

use_paper_rcparams()
panels = PanelWriter(paths.panel_dir(FIGURE))

# Row (a): three heat maps of identical geometry, so their cells line up. The top
# of each image holds the rotated reference names, the bottom the colour bar.
HEAT_SIZE = (2.13, 2.65)
HEAT_BAR = (0.12, 0.185, 0.76, 0.028)
# Row (b): the plotting box must be at least as tall as the rotated y-label is long.
MAPPER_SIZE = (7.0, 2.14)
MAPPER_AXES = (0.105, 0.215, 0.885, 0.757)

# Δ cosine: a fixed scale, not one fitted to the data -- almost every cell is near
# zero, which is the row's point, and autoscaling would turn that noise into colour.
# 0.25 is just above the largest cell (per spot -0.134 .. +0.240). The near-white
# band is stretched to ±0.1 and the ramp saturates past it, so only the references
# that really differ carry colour. Crimson = annotation ahead, teal = SpatialAIM
# ahead: hues away from the mapper colours of the dots next to the heat maps.
DELTA_LIMIT = 0.25
PALE_BAND = 0.1
DELTA_STOPS = [
    #  |delta|   annotation ahead   SpatialAIM ahead
    (0.000, "#FFFFFF", "#FFFFFF"),
    (0.050, "#FCF0F1", "#EBF7F5"),
    (PALE_BAND, "#F4C7CE", "#C7E6E1"),
    (0.175, "#B2182B", "#35978F"),
    (DELTA_LIMIT, "#67001F", "#01665E"),
]


def diverging(stops: list[tuple[float, str, str]], limit: float):
    """A colormap symmetric about 0 from stops in data units:
    (|value|, colour below zero, colour above zero), ascending from 0."""
    nodes = [(0.5 - value / (2 * limit), below) for value, below, _ in reversed(stops)]
    nodes += [(0.5 + value / (2 * limit), above) for value, _, above in stops[1:]]
    return mpl.colors.LinearSegmentedColormap.from_list("delta", nodes)


# The colour scale of each heat map: range, colormap and the colour-bar ticks.
# Δ cosine clips nothing; its inner ticks mark where the pale band ends. The
# cell-type difference (CMAP_DK, shared with figure 4) is clipped at ±DK_LIMIT.
DELTA_SCALE = {"limit": DELTA_LIMIT, "cmap": diverging(DELTA_STOPS, DELTA_LIMIT),
               "ticks": [-DELTA_LIMIT, -PALE_BAND, 0, PALE_BAND, DELTA_LIMIT],
               "labels": [f"{-DELTA_LIMIT:g}", f"{-PALE_BAND:g}", "0",
                          f"+{PALE_BAND:g}", f"+{DELTA_LIMIT:g}"]}
DK_SCALE = {"limit": DK_LIMIT, "cmap": CMAP_DK,
            "ticks": [-DK_LIMIT, 0, DK_LIMIT],
            "labels": [f"≤{-DK_LIMIT}", "0", f"≥+{DK_LIMIT}"]}


# ---------------------------------------------------------------------------
# (a) SpatialAIM against the annotation
# ---------------------------------------------------------------------------

def heat_matrix(df: pd.DataFrame, metric: str, value: str) -> pd.DataFrame:
    """The numbers of one heat map.

    Input
        df      best_per_run.csv: one row per run
        metric  the metric SpatialAIM and the annotation are compared on,
                e.g. "cossim_hard_raw_spot"
        value   what each cell shows, SpatialAIM minus annotation:
                  "delta"  the difference in ``metric``
                  "dk"     the difference in the number of cell types (K)

    Output
        a table with 10 rows (mapper x tree: nearest_centroid/ward,
        nearest_centroid/average, wann/ward, ...) and one column per reference.
        Each cell is ``value`` averaged over that reference's pairs; NaN where
        there is no pair.
    """
    pr = exp1.paired(df, metric)
    pr["dk"] = pr["k"] - pr["k_base"]
    cell = pr.groupby(["mapper", "tree", "reference"])[value].mean().unstack("reference")
    return cell.reindex(index=pd.MultiIndex.from_product([MAPPERS, TREES]),
                        columns=[r for r in REFERENCE_LABELS if r in cell.columns])


def panel_heat(matrix: pd.DataFrame, name: str, scale: dict, unit: str) -> None:
    """One heat map of row (a), its rows inside HEAT_AXES."""
    limit, cmap = scale["limit"], scale["cmap"]
    fig = plt.figure(figsize=HEAT_SIZE)
    ax = fig.add_axes(HEAT_AXES)
    ax.imshow(np.clip(matrix.to_numpy(), -limit, limit), cmap=cmap,
              vmin=-limit, vmax=limit, aspect="auto", interpolation="nearest")
    ax.xaxis.set_ticks_position("top")
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_xticklabels([REFERENCE_LABELS[r] for r in matrix.columns], rotation=90,
                       ha="center", va="bottom")
    ax.set_yticks([])
    ax.tick_params(length=0, pad=1.5)
    ax.spines[:].set_visible(False)
    for y in range(2, matrix.shape[0], 2):          # white lines between the mappers
        ax.axhline(y - 0.5, color="white", lw=1.6)

    bar = fig.add_axes(HEAT_BAR)
    bar.imshow(np.linspace(-limit, limit, 256)[None, :], cmap=cmap, aspect="auto",
               extent=(-limit, limit, 0, 1))
    bar.set_yticks([])
    bar.set_xticks(scale["ticks"])
    bar.set_xticklabels(scale["labels"])
    bar.set_xlabel(unit, labelpad=2.5)
    bar.spines[:].set_visible(False)
    bar.tick_params(length=1.5, pad=1.0)
    panels.save(fig, name)


# ---------------------------------------------------------------------------
# (b) the five mappers at their own proposal
# ---------------------------------------------------------------------------

def centred(df: pd.DataFrame, metric: str) -> pd.Series:
    """How much better or worse each mapper does than the other mappers on the same data.

    Input
        df      best_per_run.csv: one row per run; only the SpatialAIM rows
                (Leiden start, proposed K) are used
        metric  e.g. "cossim_hard_raw_spot"

    Output
        a Series indexed by (reference, mapper, tree): the mapper's ``metric``
        minus the mean of all five mappers on the same pair and tree, averaged
        over the reference's pairs. 0 = as good as the average mapper.

    Absolute values are not comparable across references (33- vs 981-gene
    panels); the offset from the other mappers on the same pair is.
    """
    runs = df[df["role"] == "spatialaim"]
    offset = runs[metric] - runs.groupby(["pair_dir", "tree"])[metric].transform("mean")
    return runs.assign(offset=offset).groupby(["reference", "mapper", "tree"])["offset"].mean()


def panel_mappers(df: pd.DataFrame) -> None:
    """Four groups (metric x tree), one box per mapper in each; a box holds the
    ten references' values."""
    groups = [(centred(df, metric), label, tree)
              for metric, label in ((METRIC_SPOT, "per spot"), (METRIC_GENE, "per gene"))
              for tree in TREES]
    fig = plt.figure(figsize=MAPPER_SIZE)
    ax = fig.add_axes(MAPPER_AXES)
    ax.axhline(0, color=PALE, lw=0.8, zorder=1)
    for g, (offsets, _, tree) in enumerate(groups):
        for m, mapper in enumerate(MAPPERS):
            colour = MAPPER_COLOURS[mapper]
            values = offsets.xs((mapper, tree), level=("mapper", "tree")).to_numpy()
            box = ax.boxplot([values], positions=[g + (m - 2) * 0.17], widths=0.15,
                             patch_artist=True, showfliers=False,
                             medianprops=dict(color=colour, lw=1.2),
                             whiskerprops=dict(color=colour, lw=0.7),
                             capprops=dict(color=colour, lw=0.7))
            box["boxes"][0].set(facecolor=colour, alpha=0.5, edgecolor=colour, linewidth=0.9)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([f"{label}\n{TREE_LABELS[tree]} linkage" for _, label, tree in groups])
    ax.set_xlim(-0.55, len(groups) - 0.45)
    ax.set_ylabel("Δ cosine similarity\nvs mean across mappers")
    bare(ax)
    panels.save(fig, "mappers")


# ---------------------------------------------------------------------------

def main() -> None:
    df = exp1.load()
    panel_heat(heat_matrix(df, METRIC_SPOT, "delta"), "delta_spot", DELTA_SCALE,
               "Δ cosine similarity per spot\n(SpatialAIM − annotation)")
    panel_heat(heat_matrix(df, METRIC_GENE, "delta"), "delta_gene", DELTA_SCALE,
               "Δ cosine similarity per gene\n(SpatialAIM − annotation)")
    panel_heat(heat_matrix(df, METRIC_SPOT, "dk"), "delta_k", DK_SCALE,
               "Δ number of cell types\n(SpatialAIM − annotation)")
    panel_mappers(df)

    panels.write_manifest()


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
