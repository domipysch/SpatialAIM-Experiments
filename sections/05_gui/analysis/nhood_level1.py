"""Squidpy nhood enrichment of the Level1 annotation mapped directly (annotation start,
nearest_centroid, K = T = 9; 6-NN, argmax states), CRC P2, as nhood_heatmap.py does for
the Leiden start at K = 16. Boxed: T cells, endothelial cells, fibroblasts (the types
SpatialAIM merges at K = 10).

Writes Data/03_GUI/Analysis/nhood_level1/: nhood_z.csv, nhood_count.csv, nhood_heatmap.{png,pdf}
"""

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import squidpy as sq
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle
from spatialaim_experiments import paths

D = paths.GUI
RUN = D / "Result/annotation_Level1"
OUT = paths.GUI_ANALYSIS / "nhood_level1"
BOXED = ["Endothelial", "Fibroblast", "T cells"]


def compute() -> tuple[pd.DataFrame, pd.DataFrame]:
    st = ad.read_h5ad(D / "ST/11_human-crc-p2.h5ad", backed="r")
    m = ad.read_h5ad(RUN / "nearest_centroid/k_009/spot_to_state_mapping_soft.h5ad")
    assert (st.obs_names == m.obs_names).all()
    # state i = start cluster i (start_cluster_to_state.csv is the identity at K = T)
    starts = ad.read_h5ad(RUN / "start_clustering.h5ad").obs
    names = (starts.drop_duplicates("start_cluster").set_index("start_cluster")
             ["start_cluster_name"].sort_index().astype(str).tolist())
    X = m.X if isinstance(m.X, np.ndarray) else m.X.toarray()
    a = ad.AnnData(obs=pd.DataFrame(index=m.obs_names))
    a.obsm["spatial"] = np.asarray(st.obsm["spatial"])
    a.obs["state"] = pd.Categorical(np.asarray(names)[X.argmax(1)], categories=names)
    print("cells per type", a.obs["state"].value_counts().reindex(names).to_dict())
    sq.gr.spatial_neighbors_knn(a, n_neighs=6)
    sq.gr.nhood_enrichment(a, cluster_key="state", n_perms=500, seed=0, n_jobs=1,
                           show_progress_bar=False)
    r = a.uns["state_nhood_enrichment"]
    z = pd.DataFrame(r["zscore"], index=names, columns=names)
    c = pd.DataFrame(r["count"], index=names, columns=names)
    z.to_csv(OUT / "nhood_z.csv")
    c.to_csv(OUT / "nhood_count.csv")
    return z, c


def heatmap(zdf: pd.DataFrame) -> plt.Figure:
    """The z-scores (symmetrised), colour clipped at ±150, the boxed types marked; no
    title (the caption says it). Also the panel of supplementary figure S11
    (SpatialAIM-Experiments/sections/05_gui/panels_s11.py)."""
    names = list(zdf.index)
    z = zdf.values
    z = (z + z.T) / 2
    n, lim = len(z), 150
    boxed = [names.index(b) for b in BOXED]
    fig, ax = plt.subplots(figsize=(7, 6.2), constrained_layout=True)
    im = ax.imshow(np.clip(z, -lim, lim), cmap="RdBu_r", norm=TwoSlopeNorm(0, -lim, lim))
    for i in range(n):
        for j in range(n):
            ax.text(j, i, f"{z[i, j]:.0f}", ha="center", va="center", fontsize=7.5,
                    color="white" if abs(z[i, j]) > 0.6 * lim else "#222")
    for i in boxed:
        for j in boxed:
            ax.add_patch(Rectangle((j - .5, i - .5), 1, 1, fill=False, ec="black", lw=1.6))
    lbl = [f"{s}*" if k in boxed else s for k, s in enumerate(names)]
    ax.set_xticks(range(n), lbl, fontsize=8, rotation=45, ha="right")
    ax.set_yticks(range(n), lbl, fontsize=8)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    fig.colorbar(im, ax=ax, shrink=0.8, label="z (colour clipped at ±150)")
    return fig


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    # the permutation test takes a while; reuse its result once written
    z_csv = OUT / "nhood_z.csv"
    fig = heatmap(pd.read_csv(z_csv, index_col=0) if z_csv.exists() else compute()[0])
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"nhood_heatmap.{ext}", dpi=200)
