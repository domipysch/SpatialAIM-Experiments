# Second row of the Results 2.2 figure: the annotation is too coarse for the panel.
#
#   (d) maps       the slice mapped from the 3-type annotation (K=3, all it can
#                  return) next to the same slice mapped from the over-clustering
#                  at the K SpatialAIM proposes
#   (e) selection  SpatialAIM's criterion scores over K for the over-clustering run
#   (f) what the three labels were summarising: each proposed state's composition
#                  in the fine annotation, which the run never saw
#
# Row letters continue from plot_invisible_split.py, so this stacks under it.
# Run from spatialaim_env (see the section README).

import json
import logging
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from spatialaim_experiments import paths
from spatialaim_experiments.style import (
    C_COHERENCE,
    C_NONN,
    C_RECON,
    C_SPATIAL,
    INK,
    PALE as C_GREY,
)

logger = logging.getLogger(__name__)

ROOT = paths.ROBUSTNESS / "granularity"
RUN_ANNOTATION = ROOT / "annotation"
RUN_OVERCLUSTER = ROOT / "overcluster"
ST_PATH = paths.DATASETS / "ST" / "01_mouse-ssp.h5ad"
OUT_DIR = paths.outputs("02_robustness")

# Which mapper's runs are read, as in the invisible-split experiment. Rebound at
# runtime by ``--mapper`` and by panels.py, so ``run_dir`` reads it rather
# than a path fixed on import.
MAPPER = "tangram"
K_ANNOTATION = 3

# Shades within a family, superficial -> deep and CGE -> MGE, so that the
# refinement of a coarse label reads as a shading of its colour. The middle of
# each ramp is the family's own colour in spatialaim_experiments.style.
C_EXC = ["#EFAE85", "#D55E00", "#7A3600"]
C_INH = ["#8FCFF0", "#0072B2", "#004E7A"]

# Dominant fine type -> the family it belongs to, and how deep it sits. The
# depth order is what makes the layer story visible in (d) and (f).
EXCITATORY = {"L2/3 IT CTX": 0, "L4/5 IT CTX": 1, "L5 IT CTX": 2, "L5 PT CTX": 2,
              "L5/6 NP CTX": 2, "L6 IT CTX": 3, "L6 CT CTX": 3, "Car3": 3, "L6b CTX": 4}
INHIBITORY = {"Vip": 0, "Sncg": 0, "Lamp5": 1, "Sst": 2, "Pvalb": 2, "Sst Chodl": 2,
              "Meis2": 1}
NON_NEURONAL = {"Endo", "Astro", "VLMC", "Oligo", "SMC-Peri", "Micro-PVM"}

# How to read a state, keyed by the fine type that dominates it -- an
# interpretation of the composition shown in (f), not an input to the run.
# The excitatory ones are named by depth, so the label does not merely
# repeat the column beneath it.
STATE_NAMES = {
    "L2/3 IT CTX": "upper layers",
    "L4/5 IT CTX": "layer 4/5",
    "L6 IT CTX": "deep layers",
    "Vip": "Vip",
    "Lamp5": "Lamp5/Sncg",
    "Sst": "Sst/Pvalb",
}
COARSE_COLOURS = {"Glutamatergic": C_EXC[1], "GABAergic": C_INH[1], "Non-Neuronal": C_NONN}

# Columns of (f): the fine types that carry the story, superficial -> deep, then
# the interneuron lineages. The six non-neuronal types are pooled -- with 1 to 10
# cells each they say nothing on their own.
FINE_ORDER = ["L2/3 IT CTX", "L4/5 IT CTX", "L5 IT CTX", "L5 PT CTX", "L5/6 NP CTX",
              "L6 IT CTX", "L6 CT CTX", "L6b CTX", "Car3",
              "Vip", "Sncg", "Lamp5", "Sst", "Pvalb", "Sst Chodl", "Meis2"]
POOLED_NONN = "non-neuronal"
# Dots below this share are dropped: single cells would otherwise speckle the panel.
MIN_DOT_SHARE = 0.02

CRITERIA = [
    ("overall", "overall", INK, 1.6),
    ("recon", "reconstruction", C_RECON, 1.0),
    ("spatial", "spatial organisation", C_SPATIAL, 1.0),
    ("coherence", "transcriptional coherence", C_COHERENCE, 1.0),
]

mpl.rcParams.update({
    "font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.0, "ytick.major.size": 2.0, "pdf.fonttype": 42,
})


def run_dir(root: Path) -> Path:
    return root / MAPPER


def cell_states(root: Path, k: int) -> tuple[np.ndarray, ad.AnnData]:
    """Reference cells labelled by the state they belong to at this cut."""
    d = run_dir(root)
    sc = ad.read_h5ad(d / "start_clustering.h5ad")
    cut = pd.read_csv(d / f"k_{k:03d}" / "start_cluster_to_state.csv")
    to_state = dict(zip(cut["start_cluster"], cut["state"]))
    return np.array([to_state[int(c)] for c in sc.obs["start_cluster"]]), sc


def spot_states(root: Path, k: int) -> np.ndarray:
    """Hard spot assignment as state ids."""
    P = np.asarray(ad.read_h5ad(run_dir(root) / f"k_{k:03d}" /
                                "spot_to_state_mapping_soft.h5ad").X)
    return np.argmax(P, axis=1)


def composition(states: np.ndarray, fine: np.ndarray) -> pd.DataFrame:
    """State x fine-type fractions, non-neuronal types pooled into one column."""
    fine = np.where(np.isin(fine, list(NON_NEURONAL)), POOLED_NONN, fine)
    table = pd.crosstab(states, fine)
    for column in FINE_ORDER + [POOLED_NONN]:
        if column not in table:
            table[column] = 0
    table = table[FINE_ORDER + [POOLED_NONN]]
    return table.div(table.sum(axis=1), axis=0)


def describe_states(fractions: pd.DataFrame) -> pd.DataFrame:
    """One row per state: family, depth rank, colour and name, from its composition."""
    rows = []
    for state, row in fractions.iterrows():
        dominant = str(row.idxmax())
        if dominant == POOLED_NONN:
            family, rank = "nonn", 0
        elif dominant in EXCITATORY:
            family, rank = "exc", EXCITATORY[dominant]
        else:
            family, rank = "inh", INHIBITORY[dominant]
        rows.append({"state": state, "dominant": dominant, "family": family, "rank": rank,
                     "name": STATE_NAMES.get(dominant, dominant)})
    out = pd.DataFrame(rows)
    order = {"exc": 0, "inh": 1, "nonn": 2}
    out = out.sort_values(["family", "rank"], key=lambda c: c.map(order) if c.name == "family" else c)
    shades = {"exc": C_EXC, "inh": C_INH, "nonn": [C_NONN]}
    colours = {}
    for family, group in out.groupby("family", sort=False):
        palette = shades[family]
        for i, state in enumerate(group["state"]):
            colours[state] = palette[min(i, len(palette) - 1)]
    out["colour"] = out["state"].map(colours)
    return out.set_index("state")


def coarse_of(column: str) -> str:
    """The coarse annotation a fine type belongs to."""
    if column == POOLED_NONN:
        return "Non-Neuronal"
    return "Glutamatergic" if column in EXCITATORY else "GABAergic"


def bare(ax):
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_map(ax, xy, colours, title):
    """Spots in random order so no state hides another."""
    bare(ax)
    order = np.random.default_rng(0).permutation(len(xy))
    ax.scatter(xy[order, 0], xy[order, 1], s=0.7, c=[colours[i] for i in order],
               linewidths=0, rasterized=True)
    ax.set_title(title, pad=2, fontsize=6.5)


def main() -> None:
    global MAPPER
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mapper", default=MAPPER,
                    choices=MAPPERS,
                    help="which mapper's runs to draw (default tangram)")
    MAPPER = ap.parse_args().mapper

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    st = ad.read_h5ad(ST_PATH, backed="r")
    xy = np.asarray(st.obsm["spatial"], dtype=float)

    best = json.loads((run_dir(RUN_OVERCLUSTER) / "k_selection.json").read_text(encoding="utf-8"))
    k_best = int(best["best"]["overall"]["k"])

    ref = ad.read_h5ad(RUN_OVERCLUSTER / "reference_scaffold.h5ad", backed="r")
    fine = ref.obs["cellTypeMinor"].astype(str).to_numpy()

    cells_over, _ = cell_states(RUN_OVERCLUSTER, k_best)
    fractions = composition(cells_over, fine)
    states = describe_states(fractions)

    # the annotation run's states carry the annotated names, so colour them by name
    _, sc_ann = cell_states(RUN_ANNOTATION, K_ANNOTATION)
    cut = pd.read_csv(run_dir(RUN_ANNOTATION) / f"k_{K_ANNOTATION:03d}" /
                      "start_cluster_to_state.csv")
    name_of = dict(zip(sc_ann.obs["start_cluster"].astype(int),
                       sc_ann.obs["start_cluster_name"].astype(str)))
    ann_colour = {int(r.state): COARSE_COLOURS.get(name_of[int(r.start_cluster)], C_GREY)
                  for r in cut.itertuples()}

    fig = plt.figure(figsize=(7.2, 3.0))

    # ---------------- (d) the slice, both starting points ----------------
    panels = [
        (RUN_ANNOTATION, K_ANNOTATION, ann_colour,
         f"from the annotation\n$K$={K_ANNOTATION}"),
        (RUN_OVERCLUSTER, k_best, states["colour"].to_dict(),
         f"from over-clustering\nproposed $K$={k_best}"),
    ]
    for i, (root, k, colours, title) in enumerate(panels):
        ax = fig.add_axes([0.030 + i * 0.150, 0.10, 0.130, 0.70])
        draw_map(ax, xy, [colours.get(s, C_GREY) for s in spot_states(root, k)], title)
    fig.text(0.008, 0.985, "d", fontsize=9, fontweight="bold", ha="left", va="top")

    # ---------------- (e) what SpatialAIM scores ----------------
    ax = fig.add_axes([0.385, 0.32, 0.200, 0.50])
    scores = pd.read_csv(run_dir(RUN_OVERCLUSTER) / "k_comparison.csv")
    ax.axvspan(k_best - 0.8, k_best + 0.8, color="0.85", alpha=0.6, lw=0)
    ax.text(k_best, 1.03, f"proposed\n$K$={k_best}", fontsize=6.5, color="0.25",
            ha="center", va="bottom")
    for key, label, colour, lw in CRITERIA:
        ax.plot(scores["k"], scores[key], lw=lw, color=colour, label=label)
        if scores[key].notna().any():
            k_star = int(scores.loc[scores[key].idxmax(), "k"])
            ax.plot([k_star], [scores[key].max()], marker="*", ms=8, color=colour,
                    mec="white", mew=0.4, zorder=5)
    ax.set_xlim(0, scores["k"].max() + 1)
    ax.set_ylim(0, 1.30)
    ax.set_xlabel("$K$"); ax.set_ylabel("criterion score")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.34), ncol=2,
              handlelength=1.2, borderpad=0.0, labelspacing=0.25, columnspacing=0.8)
    fig.text(0.335, 0.985, "e", fontsize=9, fontweight="bold", ha="left", va="top")

    # ---------------- (f) the states against the annotations ----------------
    # Rows: the seven states, unnamed -- what they are is what this panel shows.
    # Columns: the fine annotation, which the run never saw, bracketed by the
    # coarse type each fine type belongs to; a row whose dots stay inside one
    # bracket is a sub-clustering of that coarse type. Dot area is the share of
    # the row's cells, so every row sums to one; the scale is in the caption.
    ax = fig.add_axes([0.680, 0.42, 0.205, 0.38])
    columns = [c for c in fractions.columns if fractions[c].max() >= MIN_DOT_SHARE]
    rows = list(states.index)
    for x, column in enumerate(columns):
        for y, state in enumerate(rows):
            share = fractions.loc[state, column]
            if share < MIN_DOT_SHARE:
                continue
            ax.scatter([x], [y], s=3 + 80 * share, color=states.loc[state, "colour"],
                       edgecolors="0.35", linewidths=0.3, zorder=3)

    ax.set_xlim(-0.7, len(columns) - 0.3); ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels([c.replace(" CTX", "") for c in columns], rotation=90, fontsize=5.5)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([str(i + 1) for i in range(len(rows))], fontsize=6)
    for tick, state in zip(ax.get_yticklabels(), rows):
        tick.set_color(states.loc[state, "colour"])
    ax.set_xlabel("cell type in the fine annotation\n(excitatory: superficial"
                  " $\\rightarrow$ deep)", fontsize=6, labelpad=2)
    ax.set_ylabel(f"state at $K$={k_best}", fontsize=6.5, labelpad=1)

    # how each state reads, next to its row
    for y, state in enumerate(rows):
        ax.text(len(columns) - 0.1, y, states.loc[state, "name"], ha="left", va="center",
                fontsize=5.5, style="italic", color="0.3", clip_on=False)

    # bracket each run of fine types under the coarse type it belongs to
    trans = ax.get_xaxis_transform()
    groups = [coarse_of(c) for c in columns]
    start = 0
    for i in range(len(columns) + 1):
        if i == len(columns) or groups[i] != groups[start]:
            left, right = start - 0.32, i - 1 + 0.32
            middle = (left + right) / 2
            ax.plot([left, left, right, right], [1.005, 1.03, 1.03, 1.005],
                    transform=trans, color="0.55", lw=0.5, clip_on=False,
                    solid_capstyle="butt")
            if right - left > 1.5:
                ax.text(middle, 1.045, groups[start], transform=trans, ha="center",
                        va="bottom", fontsize=5.5, color="0.35")
            else:
                # a narrow group cannot carry a centred label without running into
                # its neighbour: lift it a tier and hang it off a leader line
                ax.plot([middle, middle], [1.03, 1.13], transform=trans, color="0.55",
                        lw=0.5, clip_on=False, solid_capstyle="butt")
                ax.text(middle + 0.35, 1.135, groups[start], transform=trans,
                        ha="right", va="bottom", fontsize=5.5, color="0.35")
            if i < len(columns):
                ax.axvline(i - 0.5, color="0.9", lw=0.5, zorder=1)
            start = i
    ax.grid(True, lw=0.3, color="0.92", zorder=0)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)

    fig.text(0.645, 0.985, "f", fontsize=9, fontweight="bold", ha="left", va="top")

    stem = "granularity" if MAPPER == "tangram" else f"granularity_{MAPPER}"
    for ext in ("pdf", "png"):
        fig.savefig(OUT_DIR / f"{stem}.{ext}", dpi=400)
    logger.info("Wrote %s", OUT_DIR / f"{stem}.pdf")
    logger.info("states:\n%s", states.to_string())


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
