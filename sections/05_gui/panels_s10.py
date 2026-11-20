"""
Panels of supplementary figure S10: K = 10 cell states 2 and 3 are Oliveira et al.'s
SPP1+ and SELENOP+ macrophages (Data/03_GUI, CRC P2, Leiden start, nearest centroid).

    dotplot       the GUI's custom marker dotplot of the paper's 15 macrophage genes
                  (rendered by panels.panel_dotplot, with the GUI's own function)
    registration  Visium HD tumour bins on the Xenium section after the rigid fit
    bins_spp1     their SPP1+ bins over our cell type 2 cells (50 um tumour periphery)
    bins_selenop  their SELENOP+ bins over our cell type 3 cells

and, in panels.json only (for the caption), per bin label which cell type
dominates within 24 um and the label-shuffle p-value.

The spatial check is test 2 of analysis/ (analysis/test2_markers.py, its
summary.md): the same registration (transform.json), region, counts and shuffle.
"""

from __future__ import annotations

import importlib.util
import json
import logging
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402
from spatialaim_experiments import paths  # noqa: E402
from spatialaim_experiments.style import INK, PALE, PanelWriter, use_paper_rcparams  # noqa: E402

use_paper_rcparams()

FIGURE = "S10"
ANALYSIS = Path(__file__).resolve().parent / "analysis"   # the scripts; their outputs: paths.GUI_ANALYSIS
# the GUI's tab20 colours of K = 10 cell types 2, 3 and 4 (tumour)
C2, C3, C_TUMOUR = "#ff7f0e", "#ffbb78", "#2ca02c"
C_SPP1_BIN, C_SELENOP_BIN = "#7f2704", "#08306b"   # the paper's bins: dark outlines
C_VISIUM = "#c51b7d"                                # Visium HD tumour bins
MAP_IN = (2.25, 2.55)                               # one section map, inches

panels = PanelWriter(paths.panel_dir(FIGURE), preview=paths.outputs("05_gui"))


def test2():
    """test2_markers.py as a module (its loaders and counting, not its plots)."""
    spec = importlib.util.spec_from_file_location("test2_markers", ANALYSIS / "test2_markers.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def data(t2) -> dict:
    """Test 2's inputs, with its stored registration instead of a re-fit."""
    vis = t2.load_visium()
    xen_xy, hard = t2.load_xenium()
    tr = json.loads((t2.OUT / "transform.json").read_text(encoding="utf-8"))
    reg = t2.transform(vis[["x", "y"]].to_numpy(), np.asarray(tr["centre_visium_um"]),
                       np.deg2rad(tr["angle_deg"]), tr["flip"],
                       np.asarray(tr["centre_xenium_um"]), np.asarray(tr["shift_um"]))
    mac = np.isin(hard, [t2.SPP1_STATE, t2.SELENOP_STATE])
    dist, idx = cKDTree(reg).query(xen_xy[mac], distance_upper_bound=t2.BIN_UM)
    found = np.isfinite(dist)
    periphery = (vis["Periphery"] == "50 micron").to_numpy()
    region = np.zeros(len(idx), dtype=bool)
    region[found] = periphery[idx[found]]
    sub = vis["MacrophageSubtype"].notna().to_numpy()
    return dict(vis=vis, reg=reg, xen_xy=xen_xy, hard=hard, tr=tr,
                tumour_bins=(vis["Periphery"] == "Tumor").to_numpy(),
                cells_xy=xen_xy[mac][region], state=hard[mac][region],
                bins_xy=reg[sub],
                bins_spp1=(vis.loc[sub, "MacrophageSubtype"] == "Macrophage-SPP1+").to_numpy())


def section_axes(fig, d):
    ax = fig.add_axes((0.02, 0.02, 0.96, 0.84))
    lo, hi = d["xen_xy"].min(0), d["xen_xy"].max(0)
    ax.set_xlim(lo[0] - 100, hi[0] + 100)
    ax.set_ylim(hi[1] + 100, lo[1] - 100)   # y down, as the GUI shows the section
    ax.set_aspect("equal")
    ax.axis("off")
    # 1 mm scale bar, bottom left
    x0, y0 = lo[0], hi[1] + 60
    ax.plot([x0, x0 + 1000], [y0, y0], color=INK, lw=1.2, solid_capstyle="butt")
    ax.text(x0 + 500, y0 - 70, "1 mm", ha="center", va="bottom", fontsize=7.5)
    return ax


def backdrop(ax, d, rng, colour=PALE):
    take = rng.choice(len(d["xen_xy"]), 120_000, replace=False)
    ax.scatter(*d["xen_xy"][take].T, s=0.15, c=colour, lw=0, rasterized=True)


def panel_registration(d, rng) -> None:
    fig = plt.figure(figsize=MAP_IN)
    ax = section_axes(fig, d)
    take = rng.choice(len(d["xen_xy"]), 150_000, replace=False)
    tum = d["hard"][take] == 4
    ax.scatter(*d["xen_xy"][take][~tum].T, s=0.15, c=PALE, lw=0, rasterized=True)
    ax.scatter(*d["xen_xy"][take][tum].T, s=0.15, c=C_TUMOUR, lw=0, rasterized=True)
    ax.scatter(*d["reg"][d["tumour_bins"]][::3].T, s=0.12, c=C_VISIUM, alpha=0.3, lw=0,
               rasterized=True)
    fig.text(0.02, 0.97, "Xenium tumour (cell state 4)", color=C_TUMOUR, fontsize=8, va="top")
    fig.text(0.02, 0.91, "Visium HD tumour bins, registered", color=C_VISIUM, fontsize=8,
             va="top")
    panels.save(fig, "registration")


def panel_bins(d, rng, spp1: bool) -> None:
    state, colour, edge, label, name = (
        (2, C2, C_SPP1_BIN, "SPP1+", "bins_spp1") if spp1
        else (3, C3, C_SELENOP_BIN, "SELENOP+", "bins_selenop"))
    fig = plt.figure(figsize=MAP_IN)
    ax = section_axes(fig, d)
    backdrop(ax, d, rng)
    cells = d["cells_xy"][d["state"] == state]
    bins = d["bins_xy"][d["bins_spp1"] == spp1]
    ax.scatter(*cells.T, s=1.6, c=colour, lw=0, rasterized=True)
    ax.scatter(*bins.T, s=5, marker="s", facecolors="none", edgecolors=edge, lw=0.45,
               rasterized=True)
    fig.text(0.02, 0.97, f"cell state {state} (n = {len(cells):,})", color=colour if spp1 else "#e08a2c",
             fontsize=8, va="top")
    fig.text(0.02, 0.91, f"Visium HD {label} bins (n = {len(bins):,})", color=edge,
             fontsize=8, va="top")
    panels.save(fig, name)


def bin_votes(t2, d) -> dict:
    near = t2.near_matrix(d["cells_xy"], d["bins_xy"])
    t2.RNG = np.random.default_rng(0)   # as test 2: same shuffles, same p
    votes = t2.bin_votes(near, d["state"], d["bins_spp1"])
    perm = t2.permutation_test(near, d["state"], d["bins_spp1"])
    return {"bin_votes": votes.to_dict(orient="records"), "permutation": perm.to_dict(orient="records"),
            "periphery_cells": {"state2": int((d["state"] == 2).sum()),
                                "state3": int((d["state"] == 3).sum())},
            "radius_um": t2.R_UM, "n_perm": t2.N_PERM,
            "registration": {k: d["tr"][k] for k in ("angle_deg", "flip", "shift_um", "fit")}}


def main() -> None:
    import panels as fig5

    t2 = test2()
    d = data(t2)
    fig5.panel_dotplot(panels, d["hard"])
    rng = np.random.default_rng(0)
    panel_registration(d, rng)
    panel_bins(d, rng, spp1=True)
    panel_bins(d, rng, spp1=False)
    extra = bin_votes(t2, d)
    panels.write_manifest(extra)
    for r in extra["permutation"]:
        logging.info("%s: %s %.0f%%, p = %.3g", r["bin label"], r["tested share"],
                     100 * r["observed"], r["p"])


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
