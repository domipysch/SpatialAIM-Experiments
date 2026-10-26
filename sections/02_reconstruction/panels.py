"""
Panels of the reconstruction figure (figure 2), each as its own PNG, for assembly
in Excalidraw (Excalidraw/figure_2/figure_2.excalidraw.md).

SpatialAIM is, per pair, tree and mapper, the better-scoring (score_overall) of the
Leiden start and the annotation start(s), each at its own K = best.overall
(exp1.best_start).

    a  delta_spot   Δ cosine similarity per spot, SpatialAIM - annotation (K = T)
       delta_gene   the same per gene
       delta_k      Δ number of cell types it took
                    rows: mapper x tree; columns: the ten references
    b  mappers_cell the five mappers at their own proposal, centred within each pair,
                    per cell (the spot metric)
    c  mappers_gene the same per gene, on b's y axis
    d  runtime      what each mapper costs: the K sweep from the Leiden start
    e  runtime_vs_annotation
                    a whole SpatialAIM run (Leiden scaffold + every K) against
                    mapping one annotation (its scaffold + K = T), on d's y axis

Supplementary figure 6 (figure "S6"):

    memory          peak memory of the K sweep from the Leiden start

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

FIGURE, SUPPLEMENT = "2", "S6"
LEIDEN_TAG = "leiden_res3.0"
METRIC_SPOT, METRIC_GENE = "cossim_hard_raw_spot", "cossim_hard_raw_gene"

use_paper_rcparams()
panels = PanelWriter(paths.panel_dir(FIGURE))
supplement = PanelWriter(paths.panel_dir(SUPPLEMENT))

# Row (a): three heat maps of identical geometry, so their cells line up. The top
# of each image holds the rotated reference names, the bottom the colour bar.
HEAT_SIZE = (2.13, 2.65)
HEAT_BAR = (0.12, 0.185, 0.76, 0.028)
# Row (b-e): the mapper comparison, per cell and per gene, and the two cost panels,
# one plotting band for all four (same height and vertical fractions); with three
# 0.1 in gaps as wide as row (a), 7.0 in. The plotting box must be at least as tall as the rotated y-labels are long.
# The band is 1.62 in tall (was 2.14, squashed to make room for row f-h). The
# plotting box is 1.28 in -- which is why "all K" moved from the runtime label
# into its heading in make_excalidraw.py -- and the margin below holds one line
# of x labels, since b and c say per cell / per gene in their headings.
BAND_H, BAND_BELOW, BAND_ABOVE = 1.62, 0.28, 0.06
BAND_Y = (BAND_BELOW / BAND_H, (BAND_H - BAND_BELOW - BAND_ABOVE) / BAND_H)
# The mapper comparison is split in two that share one y axis: (b) carries its
# ticks and label, (c) only the spine and ticks, set 0.05 in in. The two plotting
# boxes are equally wide, together the 4.3 in the single panel had less one gap.
MAPPER_LEFT, MAPPER_INSET = 0.735, 0.05
MAPPER_PLOT = (4.3 - 0.1 - MAPPER_LEFT - MAPPER_INSET) / 2
MAPPER_SIZES = {"mappers_cell": (MAPPER_LEFT + MAPPER_PLOT, BAND_H),
                "mappers_gene": (MAPPER_INSET + MAPPER_PLOT, BAND_H)}
COST_SIZE = (1.25, BAND_H)
COST_AXES = (0.60 / 1.25, BAND_Y[0], 0.60 / 1.25, BAND_Y[1])
# (e) shares (d)'s y axis as (c) shares (b)'s: same width as (d), its plotting box
# set in by MAPPER_INSET and nearly twice as wide, for ten boxes against d's five.
PAIR_AXES = (MAPPER_INSET / 1.25, BAND_Y[0], 1.15 / 1.25, BAND_Y[1])

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
# (b, c) the five mappers at their own proposal
# ---------------------------------------------------------------------------

def centred(df: pd.DataFrame, metric: str) -> pd.Series:
    """How much better or worse each mapper does than the other mappers on the same data.

    Input
        df      best_per_run.csv: one row per run; only the SpatialAIM rows
                (better start, K = best.overall) are used
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


def mapper_box(ax, values, position: float, width: float, mapper: str,
               outliers: bool = False, hollow: bool = False) -> None:
    """One box in the mapper's colour: median, quartiles, whiskers (1.5 x IQR);
    values beyond the whiskers as dots if ``outliers``, else not drawn. A
    ``hollow`` box is white inside, its outliers open circles."""
    colour = MAPPER_COLOURS[mapper]
    box = ax.boxplot([values], positions=[position], widths=width,
                     patch_artist=True, showfliers=outliers,
                     medianprops=dict(color=colour, lw=1.2),
                     whiskerprops=dict(color=colour, lw=0.7),
                     capprops=dict(color=colour, lw=0.7),
                     flierprops=dict(marker="o", markersize=2.5,
                                     markerfacecolor="white" if hollow else colour,
                                     markeredgecolor=colour, markeredgewidth=0.7))
    if hollow:
        box["boxes"][0].set(facecolor="white", edgecolor=colour, linewidth=0.9)
    else:
        box["boxes"][0].set(facecolor=colour, alpha=0.5, edgecolor=colour, linewidth=0.9)


def panel_mappers(df: pd.DataFrame) -> None:
    """Two panels, per cell (b) and per gene (c), on one y axis. Each has two
    groups (tree), one box per mapper in each; a box holds the ten references'
    values."""
    axes = {}
    for name, metric in (("mappers_cell", METRIC_SPOT), ("mappers_gene", METRIC_GENE)):
        w, h = MAPPER_SIZES[name]
        left = MAPPER_LEFT if name == "mappers_cell" else MAPPER_INSET
        fig = plt.figure(figsize=(w, h))
        ax = fig.add_axes((left / w, BAND_Y[0], MAPPER_PLOT / w, BAND_Y[1]))
        ax.axhline(0, color=PALE, lw=0.8, zorder=1)
        offsets = centred(df, metric)
        for g, tree in enumerate(TREES):
            for m, mapper in enumerate(MAPPERS):
                values = offsets.xs((mapper, tree), level=("mapper", "tree")).to_numpy()
                mapper_box(ax, values, g + (m - 2) * 0.17, 0.15, mapper)
        ax.set_xticks(range(len(TREES)))
        ax.set_xticklabels([f"{TREE_LABELS[tree].capitalize()} linkage" for tree in TREES])
        ax.set_xlim(-0.55, len(TREES) - 0.45)
        bare(ax)
        axes[name] = (fig, ax)
    # the shared axis: both span what either needs, and only (b) is labelled
    lows, highs = zip(*(ax.get_ylim() for _, ax in axes.values()))
    for name, (fig, ax) in axes.items():
        ax.set_ylim(min(lows), max(highs))
        if name == "mappers_cell":
            ax.set_ylabel("Δ cosine similarity\nvs mapper mean")
        else:
            ax.tick_params(axis="y", labelleft=False)
        panels.save(fig, name)


# ---------------------------------------------------------------------------
# (d, e) what each mapper costs
# ---------------------------------------------------------------------------

def cost_per_reference(res: pd.DataFrame, column: str, unit: float) -> pd.DataFrame:
    """One cost value per (reference, mapper), for the boxes of panels d and e.

    Input
        res     resources_stage_b.csv (exp1.resources("b")): one row per mapping run
        column  "wall_s" (seconds for the whole K sweep) or "peak_rss_tree_mb"
                (peak memory of the whole process tree, the aligner's own
                process included -- the same measure for every mapper)
        unit    factor to the plotted unit, e.g. 1/60 for minutes

    Output
        columns reference, mapper, value: the median over the reference's
        pairs and both trees, Leiden starts only (the full SpatialAIM pipeline).
    """
    runs = res[(res["start_tag"] == LEIDEN_TAG) & (res["status"] == "ok")]
    return (runs.groupby(["sc", "mapper"])[column].median().mul(unit)
                .rename("value").reset_index().rename(columns={"sc": "reference"}))


def whole_runs_per_reference(res: pd.DataFrame, single: pd.DataFrame) -> pd.DataFrame:
    """The two sides of panel e, one value per (reference, mapper, side), in minutes.

    Input
        res     resources_stage_b.csv (exp1.resources("b"))
        single  resources_single_k.csv (exp1.resources_single_k()): one annotation
                mapped at K = T only, per annotation unit and mapper

    Output
        columns reference, mapper, side, value:
          side "spatialaim"   a whole SpatialAIM run: Leiden scaffold build + the
                              sweep over every K (wall_with_scaffold_s), median
                              over the reference's pairs and both trees
          side "annotation"   mapping one annotation: its scaffold build + K = T,
                              median over the reference's pairs and annotation keys

    Both sides count loading, preprocessing and each K's analysis; what differs
    is the over-clustering and the K levels beyond T.

    A single-K run that missed the scaffold cache (cache_hit False: 06_mouse-retina,
    whose files on the cluster are newer than its scaffolds) rebuilt the scaffold
    inside its own timer, so its wall_s already holds the build and nothing is added.
    """
    leiden = res[(res["start_tag"] == LEIDEN_TAG) & (res["status"] == "ok")]
    single = single[single["status"] == "ok"]
    built_inside = single["cache_hit"].astype(str).str.lower() != "true"
    single = single.assign(wall_with_scaffold_s=single["wall_s"].where(
        built_inside, single["wall_s"] + single["scaffold_build_s"]))
    sides = [frame.groupby(["sc", "mapper"])["wall_with_scaffold_s"].median().div(60)
                  .rename("value").reset_index().assign(side=side)
             for side, frame in (("spatialaim", leiden), ("annotation", single))]
    return pd.concat(sides, ignore_index=True).rename(columns={"sc": "reference"})


def cost_axes(fig, rect, ylabel: str | None, log: bool):
    """The axes of one cost panel: log (runtime) or from zero (memory)."""
    ax = fig.add_axes(rect)
    if log:
        ax.set_yscale("log")
    if ylabel:
        ax.set_ylabel(ylabel)
    bare(ax)
    return ax


def mapper_columns(ax) -> None:
    """No x ticks -- the mapper colours say which box is which. After the boxes,
    which set ticks of their own."""
    ax.set_xticks([])
    ax.set_xlim(-0.6, len(MAPPERS) - 0.4)


def log_ticks(ax) -> None:
    """Plain numbers on a log axis (1, 10, 100 rather than 10^0 ...), no minor ticks."""
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_minor_locator(mpl.ticker.NullLocator())


def panel_cost(cost: pd.DataFrame, ylabel: str, log: bool):
    """Five boxes, one per mapper, each holding the ten references' values, with
    the outliers shown (a slow reference is exactly what this panel is about).
    Returns (fig, ax), unsaved."""
    fig = plt.figure(figsize=COST_SIZE)
    ax = cost_axes(fig, COST_AXES, ylabel, log)
    for m, mapper in enumerate(MAPPERS):
        mapper_box(ax, cost.loc[cost["mapper"] == mapper, "value"].to_numpy(), m, 0.62, mapper,
                   outliers=True)
    mapper_columns(ax)
    if not log:
        ax.set_ylim(bottom=0)
    return fig, ax


def panel_whole_runs(cost: pd.DataFrame):
    """Per mapper two boxes, each holding the ten references' values: filled, a
    whole SpatialAIM run; hollow, one annotation mapped. y axis labelled by (d).
    Returns (fig, ax), unsaved."""
    fig = plt.figure(figsize=COST_SIZE)
    ax = cost_axes(fig, PAIR_AXES, None, log=True)
    for m, mapper in enumerate(MAPPERS):
        for side, dx, hollow in (("spatialaim", -0.19, False), ("annotation", 0.19, True)):
            values = cost.loc[(cost["mapper"] == mapper) & (cost["side"] == side), "value"]
            mapper_box(ax, values.to_numpy(), m + dx, 0.31, mapper, outliers=True,
                       hollow=hollow)
    mapper_columns(ax)
    grey = "#7F7F7F"
    ax.legend(handles=[mpl.patches.Patch(facecolor=grey, alpha=0.5, edgecolor=grey,
                                         label="SpatialAIM"),
                       mpl.patches.Patch(facecolor="white", edgecolor=grey,
                                         label="Annotation")],
              loc="upper left", frameon=False, fontsize=6.5, handlelength=0.9,
              handleheight=0.7, borderaxespad=0.1, labelspacing=0.2, handletextpad=0.3)
    ax.tick_params(axis="y", labelleft=False)
    return fig, ax


# ---------------------------------------------------------------------------

def main() -> None:
    df = exp1.best_start(exp1.load())
    panel_heat(heat_matrix(df, METRIC_SPOT, "delta"), "delta_spot", DELTA_SCALE,
               "Δ cosine similarity per cell\n(SpatialAIM − annotation)")
    panel_heat(heat_matrix(df, METRIC_GENE, "delta"), "delta_gene", DELTA_SCALE,
               "Δ cosine similarity per gene\n(SpatialAIM − annotation)")
    panel_heat(heat_matrix(df, METRIC_SPOT, "dk"), "delta_k", DK_SCALE,
               "Δ number of cell types\n(SpatialAIM − annotation)")
    panel_mappers(df)

    runs = exp1.resources("b")
    sweep = panel_cost(cost_per_reference(runs, "wall_s", 1 / 60), "Runtime (min)", log=True)
    whole = panel_whole_runs(whole_runs_per_reference(runs, exp1.resources_single_k()))
    # the shared axis: both span what either needs, and only (d) is labelled
    lows, highs = zip(*(ax.get_ylim() for _, ax in (sweep, whole)))
    for (fig, ax), name in ((sweep, "runtime"), (whole, "runtime_vs_annotation")):
        ax.set_ylim(min(lows), max(highs))
        log_ticks(ax)
        panels.save(fig, name)
    panels.write_manifest()

    fig, _ = panel_cost(cost_per_reference(runs, "peak_rss_tree_mb", 1 / 1024),
                        "Memory used (GB)", log=False)
    supplement.save(fig, "memory")
    supplement.write_manifest()


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
