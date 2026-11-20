"""
Panels of the GUI use-case figure (figure 5, Oliveira et al. 2025 CRC P2), each as
its own PNG, for assembly in Excalidraw (Excalidraw/figure_5/figure_5.excalidraw.md).

One run: Data/03_GUI/Result, average linkage, nearest centroid, Leiden start
(+ Level1 / Level2 starts). Third draft.

Written here (the rest of figure 5 is GUI screenshots, see screenshots.py):

    genemap_<gene>  (c) the cell-type modal's reference-UMAP map of SPP1 / STAB1
    area            (d) GUI shape-features card, area row, cell types 2 / 3

and, for supplementary figure S10 (called from panels_s10.py), panel_dotplot: the
GUI's custom marker dotplot of the 15 macrophage genes of Oliveira Fig. 4c.

All are rendered by the GUI's own Plotly functions and exported with kaleido,
exactly as its Export button does.
"""

from __future__ import annotations

import logging

import anndata as ad
import numpy as np
import pandas as pd
from spatialaim_experiments import paths
from spatialaim_experiments.style import PanelWriter, png_size

FIGURE = "5"
GUI = paths.DATA / "03_GUI"
RESULT = GUI / "Result"
RUN = RESULT / "leiden" / "nearest_centroid"
SC = GUI / "scRNA" / "11_human-crc-p2.h5ad"
ST = GUI / "ST" / "11_human-crc-p2.h5ad"

K_BEST, K_SPLIT = 10, 16
SPP1, SELENOP = 2, 3                 # K = 10 states
# Oliveira Fig. 4c, in its row order
MARKERS = ["APOC1", "FBP1", "IL1RN", "FN1", "LPL", "CHI3L1", "MMP7", "SPP1",
           "MPEG1", "MS4A6A", "FGL2", "SLC40A1", "IGKC", "STAB1", "SELENOP"]

# Plotly panels: font 12-13 px, ~9 pt authored at this many px per inch.
PX_PER_IN = 125
KALEIDO_SCALE = 5

panels = PanelWriter(paths.panel_dir(FIGURE), preview=paths.outputs("05_gui"))


def hard_labels(k: int) -> np.ndarray:
    m = ad.read_h5ad(RUN / f"k_{k:03d}" / "spot_to_state_mapping_soft.h5ad")
    return np.asarray(m.X).argmax(1)


def reference() -> ad.AnnData:
    """The scaffold's obs with the Leiden start clusters, as the GUI holds it."""
    s = ad.read_h5ad(RESULT / "reference_scaffold.h5ad", backed="r")
    start = ad.read_h5ad(RESULT / "leiden" / "start_clustering.h5ad", backed="r")
    assert (start.obs_names == s.obs_names).all()
    obs = s.obs.copy()
    obs["start_cluster"] = start.obs["start_cluster"].astype(int).to_numpy()
    return ad.AnnData(obs=obs)


def save_plotly(writer: PanelWriter, fig, name: str, width_px: int, height_px: int,
                trim: bool = False) -> None:
    """Export like the GUI's Export button (kaleido PNG) and register the panel at
    PX_PER_IN, so its fonts print at the paper's size. ``trim`` crops the empty
    (transparent) border, e.g. the slack an equal-aspect plot leaves, to 4 px."""
    import io

    from PIL import Image
    from spatialaim.gui.render.common import figure_to_bytes

    fig.update_layout(width=width_px, height=height_px)
    writer.out.mkdir(parents=True, exist_ok=True)
    path = writer.out / f"{name}.png"
    path.write_bytes(figure_to_bytes(fig, "png", scale=KALEIDO_SCALE))
    if trim:
        img = Image.open(io.BytesIO(path.read_bytes()))
        pad = 4 * KALEIDO_SCALE
        l, t, r, b = img.getchannel("A").getbbox()
        box = (max(l - pad, 0), max(t - pad, 0), min(r + pad, img.width), min(b + pad, img.height))
        img.crop(box).save(path)
        width_px, height_px = (box[2] - box[0]) / KALEIDO_SCALE, (box[3] - box[1]) / KALEIDO_SCALE
    (writer.preview / f"{name}.png").write_bytes(path.read_bytes())
    writer.panels[name] = {"file": path.name,
                           "w": round(width_px / PX_PER_IN * writer.scale, 1),
                           "h": round(height_px / PX_PER_IN * writer.scale, 1)}
    logging.info("%-26s %s px (plotly x%d)", name, png_size(path), KALEIDO_SCALE)


# ---------------------------------------------------------------------------
# (d) cell area, (S10) marker dotplot -- the GUI's own figures
# ---------------------------------------------------------------------------

def panel_area(hard10: np.ndarray) -> None:
    import plotly.graph_objects as go
    from spatialaim.gui.render.common import card_layout
    from spatialaim.gui.render.shape import render_shape_features_figure

    features = pd.read_csv(RESULT / "cell_shape_features.csv", index_col=0)
    shown = np.where(np.isin(hard10, [SPP1, SELENOP]), hard10, -1)
    full = render_shape_features_figure(RUN, K_BEST, shown, features)
    # Keep the card's first row (area); the other two are eccentricity and solidity.
    fig = go.Figure([t for t in full.data if t.yaxis == "y"])
    y = full.layout.yaxis
    # y is log10(area); the GUI's 1-2-5 ticks are unevenly spaced, so the figure
    # ticks at doublings (equal steps on the log axis)
    ticks = 5 * 2.0 ** np.arange(8)    # 5 ... 640 µm²
    fig.update_yaxes(title_text=y.title.text, tickvals=np.log10(ticks),
                     ticktext=[f"{t:g}" for t in ticks])
    card_layout(fig, height=300, violinmode="overlay", violingap=0.1,
                margin=dict(l=10, r=10, t=10, b=10))
    # larger axis text than the GUI's, for the figure's print size
    names = list(dict.fromkeys(str(v) for t in fig.data for v in np.atleast_1d(t.x)))
    fig.update_xaxes(tickfont=dict(size=17), tickangle=0, tickvals=names,
                     ticktext=[as_states(n).replace("Cell state ", "Cell state<br>")
                               for n in names])
    fig.update_yaxes(tickfont=dict(size=17), title_font=dict(size=19))
    save_plotly(panels, fig, "area", 240, 300, trim=True)


def as_states(s):
    """The paper says "cell state" where the GUI says "cell type"."""
    return s.replace("Cell type", "Cell state").replace("cell type", "cell state") \
        if isinstance(s, str) else s


def relabel_states(fig) -> None:
    """as_states() on a GUI figure's category values, trace names and axis tick texts."""
    for t in fig.data:
        for a in ("x", "y"):
            v = getattr(t, a, None)
            if v is not None and len(v) and isinstance(v[0], str):
                setattr(t, a, [as_states(s) for s in v])
        if t.name:
            t.name = as_states(t.name)
    for ax in [*fig.select_xaxes(), *fig.select_yaxes()]:
        if ax.ticktext:
            ax.ticktext = [as_states(s) for s in ax.ticktext]


def panel_dotplot(writer: PanelWriter, hard10: np.ndarray) -> None:
    from spatialaim.gui.render.common import SpatialContext
    from spatialaim.gui.render.markers import render_custom_dotplot

    st_genes = ad.read_h5ad(ST, backed="r").var_names
    fig = render_custom_dotplot(reference(), RUN, K_BEST, [SPP1, SELENOP], MARKERS, str(SC),
                                spatial=SpatialContext(str(ST), hard10, None),
                                st_genes=set(st_genes))
    # Colour bar and size legend side by side instead of stacked, and dots at 0.7
    # of the GUI's size, so 22 genes fit a low panel.
    for trace in fig.data:
        trace.marker.size = np.asarray(trace.marker.size, dtype=float) * 0.7
    fig.update_traces(selector=dict(showlegend=False), marker_colorbar=dict(
        x=1.01, len=130, y=1.0, yanchor="top", title=dict(text="Mean")))
    fig.update_layout(legend=dict(x=1.19, y=1.0, yanchor="top", title=dict(text="Fraction")),
                      margin=dict(l=10, r=10, t=30, b=80))
    relabel_states(fig)
    save_plotly(writer, fig, "dotplot", 560, 240)


# ---------------------------------------------------------------------------
# (c) the cell-type modal's reference-UMAP gene maps, with larger axis text
# ---------------------------------------------------------------------------

GENE_MAPS = ["SPP1", "STAB1"]   # markers of K = 10 cell types 2 and 3
# the modal's size (CSS px, 485 x 460) less the 30 px its title took: the plot
# area keeps the modal's shape, so the UMAP is not squashed
GENE_MAP_PX = (485, 430)


def scaffold() -> ad.AnnData:
    """The reference scaffold with the Leiden start clusters, as the GUI loads it."""
    s = ad.read_h5ad(RESULT / "reference_scaffold.h5ad")
    s.obs["start_cluster"] = reference().obs["start_cluster"].to_numpy()
    return s


def panel_gene_maps() -> None:
    """The modal's "Reference — <gene> expression" map, rendered by the GUI's own
    function; here its fonts are enlarged for the figure's print size and its
    title dropped (the figure names each map above it)."""
    from spatialaim.gui.render.cell_type import render_gene_maps_figure

    sc = scaffold()
    for gene in GENE_MAPS:
        fig = render_gene_maps_figure(sc, RUN, K_BEST, gene)
        fig.update_xaxes(title_font=dict(size=20))
        fig.update_yaxes(title_font=dict(size=20))
        fig.update_traces(marker_colorbar=dict(tickfont=dict(size=20),
                                               title=dict(text="lognorm", font=dict(size=20))))
        fig.layout.annotations = ()   # (update_layout(annotations=[]) is a no-op)
        fig.update_layout(margin=dict(t=10))
        save_plotly(panels, fig, f"genemap_{gene.lower()}", *GENE_MAP_PX, trim=True)


# ---------------------------------------------------------------------------
# the numbers the text quotes: best K, Pareto-optimal Ks, the K 10 -> 16 deltas
# ---------------------------------------------------------------------------

def choose_k_numbers() -> dict:
    from spatialaim.analysis import compare_ks as scores

    table = scores.score_table(pd.read_csv(RUN / "k_comparison.csv"))
    best_all = int(scores.best_ks(table)["all"][0])
    assert best_all == K_BEST, best_all
    row = table.set_index("k")
    return {"best_k": best_all,
            "pareto_k": sorted(int(k) for k in table["k"][scores.pareto_mask(table)]),
            "delta_16_minus_10": {c.key: float(row.at[K_SPLIT, c.key] - row.at[K_BEST, c.key])
                                  for c in scores.CRITERIA}}


def main() -> None:
    hard10 = hard_labels(K_BEST)
    panel_area(hard10)
    panel_gene_maps()
    panels.write_manifest({"choose_k": choose_k_numbers()})


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
