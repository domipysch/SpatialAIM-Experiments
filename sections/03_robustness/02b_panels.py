# Panels of the Results 2.3 figure, each as its own image, for assembly in
# Excalidraw (the vault's Excalidraw/figure_3/<mapper>/figure_3_<mapper>.excalidraw.md).
#
# One mapper per invocation (``--mapper``, required); every mapper gets the same
# figure. The panels + panels.json go to Data/02_Results/robustness/plots/<mapper>/
# and are copied into the vault folder the drawing embeds them from
# (Excalidraw/figure_3/<mapper>/panels/): the drawing can only link to files
# inside the vault.
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
#                map_split_aim        the slice, SpatialAIM at the K it proposes
#                scores_split         criterion scores over K
#   too coarse   map_overcluster_aim  the slice, SpatialAIM from over-clustering at its K
#                scores_overcluster   criterion scores over K
#                dots_states_fine     each state's composition in the fine annotation
#
# The maps are cropped to the cortex: below CROP_Y the section is white matter,
# CA1 and ventricle (Codeluppi et al. 2018, Fig. 1a), which no cortical reference
# type can match, and it is a thin strip that would otherwise double the height
# of every map.
#
# Sizes are given in final print inches, at the width the figure will have in the
# paper; panels.json records what was written, so 02c_make_excalidraw.py can place
# each image at its true size.
#
# Run from spatialaim_env (see the section README).

import argparse
import json
import logging
import shutil
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_hex

from spatialaim_experiments import paths
from spatialaim_experiments.style import (
    C_COHERENCE,
    C_GLUT,
    C_RECON,
    C_SPATIAL,
    INK,
    MAPPER_LABELS as DISPLAY,
    PALE as C_GREY,
    PanelWriter,
)

logger = logging.getLogger(__name__)

mpl.rcParams.update({
    "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 9.2, "ytick.labelsize": 9.2,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.0, "ytick.major.size": 2.0,
})


# --------------------------------------------------------------------------
# the runs (01_run_robustness_experiments.ps1) and how to read them
# --------------------------------------------------------------------------

# run folders under Data/02_Results/robustness; each holds one run root per mapper
RUN_ORIGINAL = "original_annotation"           # cellType, 3 types
RUN_SPLIT = "invisible_split"                  # cellType_split, 4 types
RUN_OVERCLUSTER = "original_leiden_overcluster"
ST_PATH = paths.DATASETS / "ST" / "01_mouse-ssp.h5ad"

MAPPER: str | None = None  # set by main() from --mapper


def run_dir(run: str) -> Path:
    """One mapper's run root; its parent holds the shared reference scaffold."""
    return paths.ROBUSTNESS / run / MAPPER


def analysis_dir() -> Path:
    """Where 02a_analyse_invisible_split.py puts this mapper's numbers."""
    return paths.ROBUSTNESS / "_analysis" / MAPPER


def spot_types(run: str, k: int) -> np.ndarray:
    """Hard spot assignment as annotated cell type names."""
    d = run_dir(run)
    sc = ad.read_h5ad(d / "start_clustering.h5ad")
    cut = pd.read_csv(d / f"k_{k:03d}" / "start_cluster_to_state.csv")
    to_state = dict(zip(cut["start_cluster"], cut["state"]))
    names: dict[int, set] = {}
    for c, nm in zip(sc.obs["start_cluster"].astype(int), sc.obs["start_cluster_name"].astype(str)):
        names.setdefault(to_state[c], set()).add(nm)
    label = {s: " + ".join(sorted(v)) for s, v in names.items()}
    return np.array([label[s] for s in spot_states(run, k)], dtype=object)


def spot_states(run: str, k: int) -> np.ndarray:
    """Hard spot assignment as state ids."""
    P = np.asarray(ad.read_h5ad(run_dir(run) / f"k_{k:03d}" /
                                "spot_to_state_mapping_soft.h5ad").X)
    return np.argmax(P, axis=1)


def cell_states(run: str, k: int) -> np.ndarray:
    """Reference cells labelled by the state they belong to at this cut."""
    d = run_dir(run)
    sc = ad.read_h5ad(d / "start_clustering.h5ad")
    cut = pd.read_csv(d / f"k_{k:03d}" / "start_cluster_to_state.csv")
    to_state = dict(zip(cut["start_cluster"], cut["state"]))
    return np.array([to_state[int(c)] for c in sc.obs["start_cluster"]])


# --------------------------------------------------------------------------
# the invisible split: its cell types and their colours
# --------------------------------------------------------------------------

# The whole argument of this figure is about one annotated cell type -- the one
# that gets split, and the one an over-clustering resolves further. GABAergic and
# Non-Neuronal cells carry no part of it, so every panel that draws cells greys
# them out -- each in a grey of its own, so the panels still show that they are
# distinct types, just not ones the figure is about.
TARGET = "Glutamatergic"
A, B = f"{TARGET}_A", f"{TARGET}_B"
# The two phantom halves get hues of their own: neither is a cell type the
# reference really has, so neither may wear the Glutamatergic orange.
C_GLUT_A = "#E6C200"
C_GLUT_B = "#8C3CC8"
C_GABA = "#CFCFCF"    # light: the larger of the two background types
# mid grey: two greys of similar weight could not be told apart, and with only 25
# cells the type needs strong contrast to show at all -- but a charcoal was lost
# against the navy of the deep-layer state in (h)
C_NONN = "#7A7A7A"
TYPE_COLOURS = {TARGET: C_GLUT, A: C_GLUT_A, B: C_GLUT_B,
                "GABAergic": C_GABA, "Non-Neuronal": C_NONN}
# A state that merges two of the background types sits between their two greys.
C_GREY_MIXED = "#C2C2C2"


def colour_of(label: str) -> str:
    """A cell's or spot's colour from its type name: the two halves merged read as
    Glutamatergic, the background types each have their own grey."""
    if label in TYPE_COLOURS:
        return TYPE_COLOURS[label]
    if set(p.strip() for p in label.split("+")) == {A, B}:
        return C_GLUT
    return C_GREY_MIXED


# --------------------------------------------------------------------------
# the over-clustering: its states read against the fine annotation
# --------------------------------------------------------------------------

# The excitatory fine types, superficial -> deep: a state whose dominant fine type
# is one of them is shown, in this order.
EXCITATORY = {"L2/3 IT CTX": 0, "L4/5 IT CTX": 1, "L5 IT CTX": 2, "L5 PT CTX": 2,
              "L5/6 NP CTX": 2, "L6 IT CTX": 3, "L6 CT CTX": 3, "Car3": 3, "L6b CTX": 4}
NON_NEURONAL = {"Endo", "Astro", "VLMC", "Oligo", "SMC-Peri", "Micro-PVM"}

# How to read a state, keyed by the fine type that dominates it -- an
# interpretation of the composition panel, not an input to the run. Named by
# depth, so the label does not merely repeat the column beneath it.
STATE_NAMES = {"L2/3 IT CTX": "upper layers", "L4/5 IT CTX": "layer 4/5",
               "L6 IT CTX": "deep layers"}

# The fine annotation's types, in the composition panel's column order. A state's
# shares are taken over all of them, the six non-neuronal types pooled -- with 1
# to 10 cells each they say nothing on their own.
FINE_ORDER = ["L2/3 IT CTX", "L4/5 IT CTX", "L5 IT CTX", "L5 PT CTX", "L5/6 NP CTX",
              "L6 IT CTX", "L6 CT CTX", "L6b CTX", "Car3",
              "Vip", "Sncg", "Lamp5", "Sst", "Pvalb", "Sst Chodl", "Meis2"]
POOLED_NONN = "non-neuronal"
# Dots below this share are dropped: single cells would otherwise speckle the panel.
MIN_DOT_SHARE = 0.02

# (k_comparison.csv column, legend label, colour, line width)
CRITERIA = [
    ("overall", "overall", INK, 1.6),
    ("recon", "reconstruction", C_RECON, 1.0),
    ("spatial", "spatial", C_SPATIAL, 1.0),
    ("coherence", "coherence", C_COHERENCE, 1.0),
]


def composition(states: np.ndarray, fine: np.ndarray) -> pd.DataFrame:
    """State x fine-type fractions, non-neuronal types pooled into one column."""
    fine = np.where(np.isin(fine, list(NON_NEURONAL)), POOLED_NONN, fine)
    table = pd.crosstab(states, fine)
    for column in FINE_ORDER + [POOLED_NONN]:
        if column not in table:
            table[column] = 0
    table = table[FINE_ORDER + [POOLED_NONN]]
    return table.div(table.sum(axis=1), axis=0)


def excitatory_states(fractions: pd.DataFrame) -> pd.DataFrame:
    """The states dominated by an excitatory fine type, superficial -> deep, each
    with that type, its name and a colour along the ramp."""
    rows = []
    for state, row in fractions.iterrows():
        dominant = str(row.idxmax())
        if dominant in EXCITATORY:
            rows.append({"state": state, "dominant": dominant, "rank": EXCITATORY[dominant],
                         "name": STATE_NAMES.get(dominant, dominant)})
    out = pd.DataFrame(rows).sort_values("rank", kind="stable").set_index("state")
    out["colour"] = state_ramp(len(out))
    return out


# --------------------------------------------------------------------------
# the figure
# --------------------------------------------------------------------------

FIGURE = "3"
panels: PanelWriter | None = None  # set by main(), once the mapper is known

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
# The plotting box inside that image. The left and bottom margins hold the tick
# labels plus an axis title -- two lines of it on the left -- and the image size
# stays fixed, so the box gives up what the titles take: the score panels have to
# stay as tall as the UMAPs beside them for the first row to be flush.
SCORES_RECT = [0.34, 0.30, 0.64, 0.665]
# The last row, next to the legends. Its size follows the data rather than being
# typed: with only the excitatory states left, both the rows and the columns that
# carry them shrink, and a fixed frame would have spread three rows over the
# height of seven. Pitches and margins are print inches and reproduce the
# original frame exactly at 15 columns x 7 rows.
COL_PITCH, ROW_PITCH = 0.25, 0.155
# left and bottom also hold the axis titles; the bottom one sits under the
# longest rotated tick label ("L5/6 NP")
DOTS_MARGIN = dict(left=0.405, right=0.082, top=0.207, bottom=0.68)
DOT_SCALE = 70        # dot area per unit share, set by the row pitch above

# The figure is placed at the paper's 5.17 in text width, so the whole drawing is
# scaled by 0.71 and every point size inside a panel is authored 1/0.71 larger to
# print at its nominal size: 9.2 pt here comes out as 6.5 pt on the page.
FS_TICK = 9.2
FS_SMALL = 7.8
FS_ROW = 8.5

# the key of the cell-type colours (02c_make_excalidraw.py draws it)
TYPE_LEGEND = [("Glutamatergic", C_GLUT), ("Glut. A (Unedited)", C_GLUT_A),
               ("Glut. B (+ marker panel)", C_GLUT_B), ("GABAergic", C_GABA),
               ("Non-neuronal", C_NONN)]

# The inferred excitatory types, superficial -> deep. A ramp,
# because they are ordered by cortical depth, but a blue one: the shades of
# Glutamatergic orange they used to carry were the same hue as the annotated type
# in (a), (e), (g), and the two readings kept being confused. Blue is free --
# GABAergic is grey in every panel that draws cells now.
#
# Three anchors, interpolated to however many excitatory states the run proposes:
# Tangram's K=7 leaves exactly these three, other mappers propose a higher K and
# leave more, and the ramp has to keep reading as one.
STATE_COLOURS = ["#A6CEE3", "#2C7FB8", "#08306B"]


# The states that are not excitatory are greys in (h), one shade each: with K=7
# and only three hues, the map would otherwise look as if it had four states, not
# seven. Their greys follow the coarse greys of (a), (e): light ones around the
# GABAergic grey for a state dominated by an inhibitory fine type, the charcoal
# for a non-neuronal one.
GABA_GREYS = ("#B0B0B0", "#E2E2E2")
NONN_GREYS = (C_NONN, "#949494")


def grey_ramp(n: int, ends: tuple[str, str]) -> list[str]:
    """``n`` greys between ``ends``; the first end when ``n`` is one."""
    ramp = LinearSegmentedColormap.from_list("greys", list(ends))
    return [to_hex(ramp(0.0 if n == 1 else i / (n - 1))) for i in range(n)]


def other_state_greys(fractions: pd.DataFrame, shown) -> tuple[dict, list]:
    """A grey per state not in ``shown``, by the family of its dominant fine type,
    in state order within a family -- and the key entries for them, one per
    family present: (label, its greys)."""
    others = [s for s in fractions.index if s not in set(shown)]
    nonn = [s for s in others if fractions.loc[s].idxmax() == POOLED_NONN]
    gaba = [s for s in others if s not in nonn]
    greys, key = {}, []
    for group, ends, family in ((gaba, GABA_GREYS, "GABAergic"),
                                (nonn, NONN_GREYS, "Non-neuronal")):
        if not group:
            continue
        shades = grey_ramp(len(group), ends)
        greys.update(zip(sorted(group), shades))
        key.append((family if len(group) == 1 else f"{family} ({len(group)})", shades))
    return greys, key


def state_ramp(n: int) -> list[str]:
    """``n`` shades along the anchors above; exactly them when ``n`` is three."""
    if n <= 1:
        return STATE_COLOURS[:1]
    ramp = LinearSegmentedColormap.from_list("states", STATE_COLOURS)
    return [to_hex(ramp(i / (n - 1))) for i in range(n)]


# --------------------------------------------------------------------------
# points panels: one axes filling the whole image, so the panel is the data
# --------------------------------------------------------------------------

def points_figure(size: tuple[float, float]):
    """One frameless axes filling the whole image."""
    fig = plt.figure(figsize=size)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_aspect("equal")
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
    fig, ax = points_figure((UMAP_W, UMAP_W * dy / dx))
    draw(ax, embedding, [colour_of(x) for x in labels], 0.7)
    ax.set_xlim(x0 - margin * dx, x1 + margin * dx)
    ax.set_ylim(y0 - margin * dy, y1 + margin * dy)
    panels.save(fig, name)


def map_panel(name: str, xy, colours) -> None:
    """All four maps are the same size, so any two of them can be compared."""
    x0, x1 = float(xy[:, 0].min()), float(xy[:, 0].max())
    y1 = float(xy[:, 1].max())
    fig, ax = points_figure((MAP_W, MAP_W * (y1 - CROP_Y) / (x1 - x0)))
    draw(ax, xy, colours, 1.2)
    ax.set_xlim(x0, x1)
    ax.set_ylim(CROP_Y, y1)
    panels.save(fig, name)


# --------------------------------------------------------------------------
# scores panel: axes stay, legend does not -- both runs share one legend
# --------------------------------------------------------------------------

def scores_panel(name: str, scores: pd.DataFrame, k_best: int, dense: bool) -> None:
    fig = plt.figure(figsize=SCORES_SIZE)
    ax = fig.add_axes(SCORES_RECT)
    half = 0.8 if dense else 0.3
    ax.axvspan(k_best - half, k_best + half, color="0.85", alpha=0.6, lw=0)
    ax.text(k_best, 1.03, f"$K$={k_best}", fontsize=FS_TICK, color="0.25",
            ha="center", va="bottom")
    for key, _label, colour, lw in CRITERIA:
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
    ax.set_xlabel("$K$", fontsize=FS_TICK, labelpad=1.5)
    # Set on two lines because it is set rotated: "mapping scores" on one line is
    # 1.03 in at this size and the plotting box is 0.75 in tall, so a single line
    # would run out of both ends of the image.
    ax.set_ylabel("Mapping\nscores", fontsize=FS_TICK, labelpad=1.5,
                  linespacing=1.05)
    panels.save(fig, name)


# --------------------------------------------------------------------------
# dot panel: self-contained, its brackets and labels are column-aligned
# --------------------------------------------------------------------------

def dots_panel(name: str, fractions: pd.DataFrame,
               states: pd.DataFrame) -> list[float]:
    """``fractions`` is already restricted to the states being shown, so the
    columns follow: a fine type with no mass in any of them drops out of the
    panel entirely rather than sitting there as an empty column."""
    columns = [c for c in fractions.columns if fractions[c].max() >= MIN_DOT_SHARE]
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
            if share < MIN_DOT_SHARE:
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
    ax.set_xlabel("Fine cell type (reference)", fontsize=FS_TICK, labelpad=3)
    ax.set_ylabel("Found type", fontsize=FS_TICK, labelpad=4)

    # one bracket over all columns: they are the coarse type the figure is about
    trans = ax.get_xaxis_transform()
    left, right = -0.32, len(columns) - 1 + 0.32
    ax.plot([left, left, right, right], [1.015, 1.05, 1.05, 1.015],
            transform=trans, color="0.55", lw=0.5, clip_on=False, solid_capstyle="butt")
    ax.text((left + right) / 2, 1.06, TARGET, transform=trans, ha="center",
            va="bottom", fontsize=FS_SMALL, color="0.35")
    ax.grid(True, lw=0.3, color="0.92", zorder=0)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0, pad=1.5)
    panels.save(fig, name)
    # where each row sits in the panel, as a fraction of its height from the top
    bottom, height = axes[1], axes[3]
    n = len(rows)
    return [1 - bottom - height * (n - 0.4 - i) / (n + 0.2) for i in range(n)]


def main() -> None:
    global panels, MAPPER
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mapper", required=True, choices=sorted(DISPLAY),
                    help="which mapper's runs to draw")
    mapper = MAPPER = ap.parse_args().mapper
    out = paths.ROBUSTNESS / "plots" / mapper
    panels = PanelWriter(out)

    st = ad.read_h5ad(ST_PATH, backed="r")
    xy = np.asarray(st.obsm["spatial"], dtype=float)

    # ---------------- the reference and its two annotations ----------------
    scaffold = ad.read_h5ad(run_dir(RUN_SPLIT).parent / "reference_scaffold.h5ad", backed="r")
    emb = np.asarray(scaffold.obsm["X_umap_shared_genes"])
    umap_panel("umap_celltype", emb, scaffold.obs["cellType"].astype(str).to_numpy())
    umap_panel("umap_split", emb, scaffold.obs["cellType_split"].astype(str).to_numpy())

    res = json.loads((analysis_dir() / "results.json").read_text(encoding="utf-8"))
    # (g) is whatever level this run proposes -- K=3, the re-merge, for every
    # mapper so far, but the panel follows the proposal rather than asserting it
    k_split = int(res["k_proposed"]["overall"])
    for name, run, k in [("map_celltype", RUN_ORIGINAL, 3),
                         ("map_split", RUN_SPLIT, 4),
                         ("map_split_aim", RUN_SPLIT, k_split)]:
        map_panel(name, xy, [colour_of(t) for t in spot_types(run, k)])

    scores_panel("scores_split", pd.DataFrame(res["scores"]), k_split, dense=False)

    # ---------------- too coarse: over-clustering as the starting point ------
    selection = json.loads((run_dir(RUN_OVERCLUSTER) / "k_selection.json")
                           .read_text(encoding="utf-8"))
    k_best = int(selection["best"]["overall"]["k"])

    ref = ad.read_h5ad(run_dir(RUN_OVERCLUSTER).parent / "reference_scaffold.h5ad", backed="r")
    fine = ref.obs["cellTypeMinor"].astype(str).to_numpy()
    fractions = composition(cell_states(RUN_OVERCLUSTER, k_best), fine)
    # Only the excitatory states are shown: the other families are grey on the
    # slice and absent from the composition panel, so the figure asks one
    # question throughout.
    states = excitatory_states(fractions)

    greys, other_key = other_state_greys(fractions, states.index)
    colours = {**greys, **states["colour"].to_dict()}
    map_panel("map_overcluster_aim", xy,
              [colours.get(s, C_GREY)
               for s in spot_states(RUN_OVERCLUSTER, k_best)])
    scores_panel("scores_overcluster",
                 pd.read_csv(run_dir(RUN_OVERCLUSTER) / "k_comparison.csv"),
                 k_best, dense=True)
    # Columns too, not just rows: state 1 carries a little Lamp5, which would
    # otherwise leave one GABAergic column standing alone in a panel that is
    # about the excitatory types.
    excitatory = [c for c in fractions.columns if c in EXCITATORY]
    row_y = dots_panel("dots_states_fine",
                       fractions.loc[states.index, excitatory], states)

    # ---------------- what the legends need to say ----------------
    panels.write_manifest({
        "mapper": mapper,
        "mapper_name": DISPLAY[mapper],
        "k_split": k_split,
        "k_overcluster": k_best,
        "legends": {
            "annotation": TYPE_LEGEND,
            "criteria": [[label, colour, lw] for _key, label, colour, lw in CRITERIA],
            "states": [{"n": i + 1, "colour": states.loc[s, "colour"],
                        "name": states.loc[s, "name"],
                        "dominant": states.loc[s, "dominant"], "y": row_y[i]}
                       for i, s in enumerate(states.index)],
            # the grey states of (h), keyed under the excitatory ones
            "other_states": [[label, shades] for label, shades in other_key],
        },
    })
    logger.info("states:\n%s", states.to_string())

    # the copy the Excalidraw drawing embeds
    vault_dir = paths.panel_dir(FIGURE, mapper)
    shutil.copytree(out, vault_dir, dirs_exist_ok=True)
    logger.info("copied to %s", vault_dir)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
