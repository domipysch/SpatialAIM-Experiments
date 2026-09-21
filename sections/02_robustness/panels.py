# Panels of the Results 2.2 figure, each as its own image, for assembly in
# Excalidraw (C:\Users\zi69hebi\Dev\00_Obsidian\Excalidraw\figure_2_2.excalidraw.md).
#
# One mapper per invocation (``--mapper``, default Tangram, which is the figure
# the paper prints). Every other mapper writes its panels to
# Excalidraw\figure_2_2_<mapper>\ and is assembled into its own drawing, so the
# same figure exists once per mapper.
#
# The figure carries almost no text: what is not anchored to data (letters,
# legends) is drawn in Excalidraw, and everything else belongs in the caption.
# What stays in these images is the points, the axes, the grey band marking the
# proposed K, and the composition panel's brackets and tick labels.
#
#   too fine     umap_celltype        reference in shared-gene space, coarse annotation
#                umap_split           the same embedding, the split annotation
#                map_celltype         the slice, Tangram on the three annotated types
#                map_split            the slice, Tangram on the split annotation
#                map_split_aim        the slice, AIM at the K it proposes
#                scores_split         criterion scores over K
#   too coarse   map_overcluster_aim  the slice, AIM from over-clustering at its K
#                scores_overcluster   criterion scores over K
#                dots_states_fine     each state's composition in the fine annotation
#
# The maps are cropped to the cortex: below CROP_Y the section is white matter,
# CA1 and ventricle (Codeluppi et al. 2018, Fig. 1a), which no cortical reference
# type can match, and it is a thin strip that would otherwise double the height
# of every map.
#
# Sizes are given in final print inches, at the width the figure will have in the
# paper; panels.json records what was written, so make_excalidraw.py can place
# each image at its true size.
#
# Run from aim_env (see the section README).

import argparse
import json
import logging

import anndata as ad
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_hex

import plot_granularity as gran               # sibling scripts: the data access
import plot_invisible_split as split          # and palette of each half

from aim_experiments import paths
from aim_experiments.style import PanelWriter

logger = logging.getLogger(__name__)

# both plot modules set their own rcParams on import; this figure needs larger
plt.rcParams.update({"xtick.labelsize": 9.2, "ytick.labelsize": 9.2})

FIGURE = "2_2"
# rebound by main() once the mapper is known; Tangram keeps the unsuffixed folder
panels = PanelWriter(paths.panel_dir(FIGURE))

# the cortex ends here: below it the section narrows into white matter, CA1 and
# the ventricle
CROP_Y = 23000.0

# Four even columns, as wide as the row can afford: 4 x 470 units plus three 84
# unit gaps fill the figure's width, and nothing lives between the columns any
# more. Every panel of a row is the same width, and the score panels are as tall
# as the UMAPs so the first row is flush.
COL_W = 470 / 300
UMAP_W = COL_W        # the embedding is wider than tall, so the panel is too
MAP_W = COL_W
SCORES_SIZE = (COL_W, 337 / 300)
# The last row, next to the legends. Its size follows the data rather than being
# typed: with only the excitatory states left, both the rows and the columns that
# carry them shrink, and a fixed frame would have spread three rows over the
# height of seven. Pitches and margins are print inches and reproduce the
# original frame exactly at 15 columns x 7 rows.
COL_PITCH, ROW_PITCH = 0.25, 0.155
DOTS_MARGIN = dict(left=0.245, right=0.082, top=0.207, bottom=0.775)
DOT_SCALE = 70        # dot area per unit share, set by the row pitch above

# The figure is placed at the paper's 5.17 in text width, so the whole drawing is
# scaled by 0.71 and every point size inside a panel is authored 1/0.71 larger to
# print at its nominal size: 9.2 pt here comes out as 6.5 pt on the page.
FS_TICK = 9.2
FS_SMALL = 7.8
FS_ROW = 8.5

# criterion names short enough for a legend column
CRITERION_LABELS = {"recon": "reconstruction", "spatial": "spatial",
                    "coherence": "coherence", "overall": "overall"}
COARSE_LEGEND = [("Glutamatergic", split.C_GLUT), ("GABAergic", split.C_GABA),
                 ("Non-Neuronal", split.C_NONN)]
SPLIT_LEGEND = [("Glutamatergic A", split.C_GLUT_A), ("Glutamatergic B", split.C_GLUT_B)]

# The whole argument of this figure is about one annotated cell type -- the one
# that gets split, and the one an over-clustering resolves further. GABAergic and
# Non-Neuronal cells carry no part of it, so in every panel that draws cells they
# are greyed out rather than coloured: (a), (b), (e), (f), (g) by label, and (h)
# by the family of each AIM state's dominant fine type. The composition panel (i)
# keeps all three families in colour -- it is the one place the other two are the
# point.
GREYED_TYPES = {"GABAergic", "Non-Neuronal"}
FOCUS_FAMILIES = {"exc"}

# The inferred types that survive that focus, superficial -> deep. A ramp,
# because they are ordered by cortical depth, but a blue one: the shades of
# Glutamatergic orange they used to carry were the same hue as the annotated type
# in (a), (e), (g), and the two readings kept being confused. Blue is free --
# GABAergic is grey in every panel that draws cells now.
#
# Three anchors, interpolated to however many excitatory states the run proposes:
# Tangram's K=7 leaves exactly these three, other mappers propose a higher K and
# leave more, and the ramp has to keep reading as one.
STATE_COLOURS = ["#A6CEE3", "#2C7FB8", "#08306B"]


def state_ramp(n: int) -> list[str]:
    """``n`` shades along the anchors above; exactly them when ``n`` is three."""
    if n <= 1:
        return STATE_COLOURS[:1]
    ramp = LinearSegmentedColormap.from_list("states", STATE_COLOURS)
    return [to_hex(ramp(i / (n - 1))) for i in range(n)]


def focus_colour(label: str) -> str:
    """``split.colour_of``, with the two families this figure is not about grey."""
    return split.C_GREY if label in GREYED_TYPES else split.colour_of(label)

def save(fig, name: str) -> None:
    panels.save(fig, name)


# --------------------------------------------------------------------------
# points panels: one axes filling the whole image, so the panel is the data
# --------------------------------------------------------------------------

def points_figure(size: tuple[float, float], frame: bool):
    fig = plt.figure(figsize=size)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_aspect("equal")
    for spine in ax.spines.values():
        spine.set_visible(frame)
        spine.set_linewidth(0.5)
        spine.set_edgecolor("0.75")
    return fig, ax


def draw(ax, coords, colours, size: float) -> None:
    """Points in random order so no group hides another."""
    order = np.random.default_rng(0).permutation(len(coords))
    ax.scatter(coords[order, 0], coords[order, 1], s=size,
               c=[colours[i] for i in order], linewidths=0, rasterized=True)


def umap_panel(name: str, embedding, labels) -> None:
    """Cropped to the embedding's own extent: a square frame would be mostly
    empty, and the panel has a column's width to fill."""
    margin = 0.03
    (x0, y0), (x1, y1) = embedding.min(axis=0), embedding.max(axis=0)
    dx, dy = (x1 - x0) * (1 + 2 * margin), (y1 - y0) * (1 + 2 * margin)
    fig, ax = points_figure((UMAP_W, UMAP_W * dy / dx), frame=True)
    draw(ax, embedding, [focus_colour(x) for x in labels], 0.7)
    ax.set_xlim(x0 - margin * dx, x1 + margin * dx)
    ax.set_ylim(y0 - margin * dy, y1 + margin * dy)
    save(fig, name)


def map_panel(name: str, xy, colours) -> None:
    """All four maps are the same size, so any two of them can be compared."""
    x0, x1 = float(xy[:, 0].min()), float(xy[:, 0].max())
    y1 = float(xy[:, 1].max())
    fig, ax = points_figure((MAP_W, MAP_W * (y1 - CROP_Y) / (x1 - x0)), frame=False)
    draw(ax, xy, colours, 1.2)
    ax.set_xlim(x0, x1)
    ax.set_ylim(CROP_Y, y1)
    save(fig, name)


# --------------------------------------------------------------------------
# scores panel: axes stay, legend does not -- both runs share one legend
# --------------------------------------------------------------------------

def scores_panel(name: str, scores: pd.DataFrame, k_best: int, dense: bool) -> None:
    fig = plt.figure(figsize=SCORES_SIZE)
    ax = fig.add_axes([0.30, 0.235, 0.68, 0.735])
    half = 0.8 if dense else 0.3
    ax.axvspan(k_best - half, k_best + half, color="0.85", alpha=0.6, lw=0)
    ax.text(k_best, 1.03, f"$K$={k_best}", fontsize=FS_TICK, color="0.25",
            ha="center", va="bottom")
    for key, _label, colour, lw in gran.CRITERIA:
        if not scores[key].notna().any():
            continue
        if dense:
            ax.plot(scores["k"], scores[key], lw=lw, color=colour)
        else:
            ax.plot(scores["k"], scores[key], marker="o", ms=2.6, lw=lw, color=colour)
        k_star = int(scores.loc[scores[key].idxmax(), "k"])
        ax.plot([k_star], [scores[key].max()], marker="*", ms=7, color=colour,
                mec="white", mew=0.4, zorder=5)
    if dense:
        ax.set_xlim(0, int(scores["k"].max()) + 1)
    else:
        ax.set_xticks(sorted(scores["k"].astype(int)))
    ax.set_yticks([0.0, 0.5, 1.0])
    ax.set_ylim(0, 1.25)
    save(fig, name)


# --------------------------------------------------------------------------
# dot panel: self-contained, its brackets and labels are column-aligned
# --------------------------------------------------------------------------

def dots_panel(name: str, fractions: pd.DataFrame,
               states: pd.DataFrame) -> list[float]:
    """``fractions`` is already restricted to the states being shown, so the
    columns follow: a fine type with no mass in any of them drops out of the
    panel entirely rather than sitting there as an empty column."""
    columns = [c for c in fractions.columns if fractions[c].max() >= gran.MIN_DOT_SHARE]
    rows = list(states.index)
    m = DOTS_MARGIN
    size = (m["left"] + len(columns) * COL_PITCH + m["right"],
            m["top"] + len(rows) * ROW_PITCH + m["bottom"])
    axes = [m["left"] / size[0], m["bottom"] / size[1],
            len(columns) * COL_PITCH / size[0], len(rows) * ROW_PITCH / size[1]]
    fig = plt.figure(figsize=size)
    ax = fig.add_axes(axes)
    for x, column in enumerate(columns):
        for y, state in enumerate(rows):
            share = fractions.loc[state, column]
            if share < gran.MIN_DOT_SHARE:
                continue
            ax.scatter([x], [y], s=2 + DOT_SCALE * share,
                       color=states.loc[state, "colour"],
                       edgecolors="0.35", linewidths=0.3, zorder=3)

    ax.set_xlim(-0.7, len(columns) - 0.3); ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels([c.replace(" CTX", "") for c in columns], rotation=90, fontsize=FS_SMALL)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([str(i + 1) for i in range(len(rows))], fontsize=FS_ROW,
                       color="#1e1e1e")

    trans = ax.get_xaxis_transform()
    groups = [gran.coarse_of(c) for c in columns]
    start = 0
    for i in range(len(columns) + 1):
        if i == len(columns) or groups[i] != groups[start]:
            left, right = start - 0.32, i - 1 + 0.32
            middle = (left + right) / 2
            ax.plot([left, left, right, right], [1.015, 1.05, 1.05, 1.015],
                    transform=trans, color="0.55", lw=0.5, clip_on=False,
                    solid_capstyle="butt")
            if right - left > 1.5:
                ax.text(middle, 1.06, groups[start], transform=trans, ha="center",
                        va="bottom", fontsize=FS_SMALL, color="0.35")
            else:
                # a narrow group cannot carry a centred label without running into
                # its neighbour: lift it a tier and hang it off a leader line
                ax.plot([middle, middle], [1.05, 1.105], transform=trans, color="0.55",
                        lw=0.5, clip_on=False, solid_capstyle="butt")
                ax.text(middle + 0.35, 1.11, groups[start], transform=trans,
                        ha="right", va="bottom", fontsize=FS_SMALL, color="0.35")
            if i < len(columns):
                ax.axvline(i - 0.5, color="0.9", lw=0.5, zorder=1)
            start = i
    ax.grid(True, lw=0.3, color="0.92", zorder=0)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0, pad=1.5)
    save(fig, name)
    # where each row sits in the panel, as a fraction of its height from the top
    bottom, height = axes[1], axes[3]
    n = len(rows)
    return [1 - bottom - height * (n - 0.4 - i) / (n + 0.2) for i in range(n)]


def main() -> None:
    global panels
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mapper", default="tangram", choices=sorted(split.DISPLAY),
                    help="which mapper's runs to draw (default tangram)")
    mapper = ap.parse_args().mapper
    # both plot modules resolve their run directories through this
    gran.MAPPER = split.MAPPER = mapper
    panels = PanelWriter(paths.panel_dir(FIGURE,
                                         None if mapper == "tangram" else mapper))

    st = ad.read_h5ad(gran.ST_PATH, backed="r")
    xy = np.asarray(st.obsm["spatial"], dtype=float)

    # ---------------- the reference and its two annotations ----------------
    scaffold = ad.read_h5ad(split.ROOT / "split" / "reference_scaffold.h5ad", backed="r")
    emb = np.asarray(scaffold.obsm["X_umap_shared_genes"])
    umap_panel("umap_celltype", emb, scaffold.obs["cellType"].astype(str).to_numpy())
    umap_panel("umap_split", emb, scaffold.obs["cellType_split"].astype(str).to_numpy())

    res = json.loads((split.analysis_dir() / "results.json").read_text(encoding="utf-8"))
    # (g) is whatever level this run proposes -- K=3, the re-merge, for every
    # mapper so far, but the panel follows the proposal rather than asserting it
    k_split = int(res["selection"]["split"]["best"]["overall"])
    for name, reference, k in [("map_celltype", "original", 3),
                               ("map_split", "split", 4),
                               ("map_split_aim", "split", k_split)]:
        types = split.spot_types(reference, k)
        map_panel(name, xy, [focus_colour(t) for t in types])

    scores_panel("scores_split", pd.DataFrame(res["scores"]["split"]), k_split, dense=False)

    # ---------------- too coarse: over-clustering as the starting point ------
    selection = json.loads((gran.run_dir(gran.RUN_OVERCLUSTER) / "k_selection.json")
                           .read_text(encoding="utf-8"))
    k_best = int(selection["best"]["overall"]["k"])

    ref = ad.read_h5ad(gran.RUN_OVERCLUSTER / "reference_scaffold.h5ad", backed="r")
    fine = ref.obs["cellTypeMinor"].astype(str).to_numpy()
    cells, _ = gran.cell_states(gran.RUN_OVERCLUSTER, k_best)
    fractions = gran.composition(cells, fine)
    states = gran.describe_states(fractions)

    # Everything below this line is about the excitatory states alone: the other
    # families are grey on the slice and absent from the composition panel, so
    # the figure asks one question throughout.
    states = states[states["family"].isin(FOCUS_FAMILIES)].copy()
    states["colour"] = state_ramp(len(states))

    colours = states["colour"].to_dict()
    map_panel("map_overcluster_aim", xy,
              [colours.get(s, gran.C_GREY)
               for s in gran.spot_states(gran.RUN_OVERCLUSTER, k_best)])
    scores_panel("scores_overcluster",
                 pd.read_csv(gran.run_dir(gran.RUN_OVERCLUSTER) / "k_comparison.csv"),
                 k_best, dense=True)
    # Columns too, not just rows: state 1 carries a little Lamp5, which would
    # otherwise leave one GABAergic column and its bracket standing alone in a
    # panel that is about the excitatory types.
    excitatory = [c for c in fractions.columns if gran.coarse_of(c) == "Glutamatergic"]
    row_y = dots_panel("dots_states_fine",
                       fractions.loc[states.index, excitatory], states)

    # ---------------- what the legends need to say ----------------
    panels.write_manifest({
        "mapper": mapper,
        "mapper_name": split.DISPLAY[mapper],
        "k_split": k_split,
        "k_overcluster": k_best,
        "legends": {
            "annotation": COARSE_LEGEND[:1] + SPLIT_LEGEND + COARSE_LEGEND[1:],
            "criteria": [[CRITERION_LABELS[key], colour, lw]
                         for key, _label, colour, lw in gran.CRITERIA],
            "states": [{"n": i + 1, "colour": states.loc[s, "colour"],
                        "name": states.loc[s, "name"],
                        "dominant": states.loc[s, "dominant"], "y": row_y[i]}
                       for i, s in enumerate(states.index)],
        },
    })
    logger.info("states:\n%s", states.to_string())


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
