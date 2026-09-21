# Figure for the invisible-split half of Results 2.2, for one mapper
# (``--mapper``, default Tangram).
#
#   (a) premise   the reference in its shared-gene embedding, coloured by the
#                 original annotation and by the annotation carrying the split --
#                 same embedding twice, because the panel cannot tell the halves
#                 apart
#   (b) maps      the mapper on the original annotation, the mapper on the
#                 annotation with the split, and SpatialAIM's output at the K it proposes
#   (c) selection SpatialAIM's criterion scores over K, with the K each one picks
#
# Run from spatialaim_env (see the section README).

import argparse
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
    C_GABA,
    C_GLUT,
    C_NONN,
    C_RECON,
    C_SPATIAL,
    INK,
    MAPPER_LABELS as DISPLAY,
    PALE as C_GREY,
)

logger = logging.getLogger(__name__)

ROOT = paths.ROBUSTNESS / "invisible_split"
ST_PATH = paths.DATASETS / "ST" / "01_mouse-ssp.h5ad"
OUT_DIR = paths.outputs("02_robustness")

# Which spot-to-state mapper's runs are read. Rebound by the callers that draw
# the figure for another mapper (panels.py), so everything below reaches it
# through the functions rather than through a constant computed on import.
MAPPER = "tangram"
TARGET = "Glutamatergic"
A, B = f"{TARGET}_A", f"{TARGET}_B"

# The two phantom halves get hues of their own: neither is a cell type the
# reference really has, so neither may wear the Glutamatergic orange.
C_GLUT_A = "#E6C200"
C_GLUT_B = "#8C3CC8"

TYPE_COLOURS = {
    "GABAergic": C_GABA,
    TARGET: C_GLUT,
    A: C_GLUT_A,
    B: C_GLUT_B,
    "Non-Neuronal": C_NONN,
}
PRETTY = {A: "Glutamatergic A", B: "Glutamatergic B"}
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


def colour_of(label: str) -> str:
    if label in TYPE_COLOURS:
        return TYPE_COLOURS[label]
    if set(p.strip() for p in label.split("+")) == {A, B}:
        return C_GLUT
    return C_GREY


def run_dir(reference: str) -> Path:
    return ROOT / reference / MAPPER


def analysis_dir() -> Path:
    """Where analyse_invisible_split.py puts this mapper's numbers."""
    return ROOT / "_analysis" / MAPPER


def spot_types(reference: str, k: int) -> np.ndarray:
    """Hard spot assignment as annotated cell type names."""
    d = run_dir(reference)
    sc = ad.read_h5ad(d / "start_clustering.h5ad")
    cut = pd.read_csv(d / f"k_{k:03d}" / "start_cluster_to_state.csv")
    to_state = dict(zip(cut["start_cluster"], cut["state"]))
    names: dict[int, set] = {}
    for c, nm in zip(sc.obs["start_cluster"].astype(int), sc.obs["start_cluster_name"].astype(str)):
        names.setdefault(to_state[c], set()).add(nm)
    label = {s: " + ".join(sorted(v)) for s, v in names.items()}
    P = np.asarray(ad.read_h5ad(d / f"k_{k:03d}" / "spot_to_state_mapping_soft.h5ad").X)
    return np.array([label[s] for s in np.argmax(P, axis=1)], dtype=object)


def bare(ax, frame: bool = False):
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(frame)
        s.set_linewidth(0.5)


def draw_points(ax, coords, labels, title, size=0.7, frame=False):
    """Points coloured by cell type, drawn in random order so none hides another."""
    bare(ax, frame=frame)
    order = np.random.default_rng(0).permutation(len(coords))
    ax.scatter(coords[order, 0], coords[order, 1], s=size,
               c=[colour_of(labels[i]) for i in order], linewidths=0, rasterized=True)
    if title:
        ax.set_title(title, pad=2)


def main() -> None:
    global MAPPER
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mapper", default=MAPPER, choices=sorted(DISPLAY),
                    help="which mapper's runs to draw (default tangram)")
    MAPPER = ap.parse_args().mapper
    name = DISPLAY[MAPPER]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    res = json.loads((analysis_dir() / "results.json").read_text(encoding="utf-8"))
    st = ad.read_h5ad(ST_PATH, backed="r")
    xy = np.asarray(st.obsm["spatial"], dtype=float)

    fig = plt.figure(figsize=(7.2, 3.0))

    # ---------------- (a) premise ----------------
    scaffold = ad.read_h5ad(run_dir("split").parent / "reference_scaffold.h5ad", backed="r")
    emb = np.asarray(scaffold.obsm["X_umap_shared_genes"])
    # the same embedding twice: only the labels differ, because the panel cannot
    # tell the two halves apart
    for i, (column, title) in enumerate([
        ("cellType", "Original annotation"),
        ("cellType_split", "Glutamatergic type\nrandomly splitted"),
    ]):
        ax = fig.add_axes([0.050, 0.600 - i * 0.385, 0.135, 0.30])
        draw_points(ax, emb, scaffold.obs[column].astype(str).to_numpy(), title,
                    size=0.5, frame=True)
        for s in ax.spines.values():
            s.set_edgecolor("0.75")
    # that both panels are the same embedding is said in the caption, not here
    fig.text(0.008, 0.985, "a", fontsize=9, fontweight="bold", ha="left", va="top")

    handles = [
        plt.Line2D([], [], marker="o", ls="", ms=3, mfc=TYPE_COLOURS[t], mec="none",
                   label=PRETTY.get(t, t))
        for t in [TARGET, A, B, "GABAergic", "Non-Neuronal"]
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.022, 0.175), frameon=False,
               ncol=2, handletextpad=0.15, columnspacing=0.5, borderpad=0.0, labelspacing=0.25)

    # ---------------- (b) the slice ----------------
    k_split = int(res["selection"]["split"]["best"]["overall"])
    panels = [
        ("original", 3, f"{name} mapping,\noriginal annotation"),
        ("split", 4, f"{name} mapping,\nsplit annotation"),
        ("split", k_split, f"SpatialAIM with {name},\nsplit annotation,\n"
                           f"proposed solution $K$={k_split}"),
    ]
    for i, (reference, k, title) in enumerate(panels):
        ax = fig.add_axes([0.285 + i * 0.145, 0.07, 0.135, 0.76])
        draw_points(ax, xy, spot_types(reference, k), title)
        ax.title.set_fontsize(6.5)
    fig.text(0.283, 0.985, "b", fontsize=9, fontweight="bold", ha="left", va="top")

    # ---------------- (c) what SpatialAIM scores ----------------
    ax = fig.add_axes([0.800, 0.40, 0.175, 0.44])
    split = pd.DataFrame(res["scores"]["split"])
    best_overall = res["selection"]["split"]["best"]["overall"]
    ax.axvspan(best_overall - 0.3, best_overall + 0.3, color="0.85", alpha=0.6, lw=0)
    ax.text(best_overall, 1.03, f"proposed\n$K$={best_overall}", fontsize=6.5, color="0.25",
            ha="center", va="bottom")
    for key, label, colour, lw in CRITERIA:
        ax.plot(split["k"], split[key], marker="o", ms=2.6, lw=lw, color=colour, label=label)
        if split[key].notna().any():
            best = int(split.loc[split[key].idxmax(), "k"])
            ax.plot([best], [split[key].max()], marker="*", ms=8, color=colour,
                    mec="white", mew=0.4, zorder=5)
    ax.set_xticks(sorted(split["k"].astype(int)))
    ax.set_xlabel("$K$"); ax.set_ylabel("criterion score")
    ax.set_ylim(0, 1.30)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=1,
              handlelength=1.2, borderpad=0.0, labelspacing=0.25)
    fig.text(0.745, 0.985, "c", fontsize=9, fontweight="bold", ha="left", va="top")

    stem = "invisible_split" if MAPPER == "tangram" else f"invisible_split_{MAPPER}"
    for ext in ("pdf", "png"):
        fig.savefig(OUT_DIR / f"{stem}.{ext}", dpi=400)
    logger.info("Wrote %s", OUT_DIR / f"{stem}.pdf")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
