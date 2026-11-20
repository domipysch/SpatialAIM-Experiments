"""Heatmaps of squidpy nhood enrichment, CRC P2, Leiden start, nearest_centroid, K=16 (6-NN, argmax states)."""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle
from spatialaim_experiments import paths
D = paths.GUI_ANALYSIS / "nhood_k016"
G = [6, 8, 9, 10]
z = pd.read_csv(D / "nhood_z.csv", index_col=0).values
c = pd.read_csv(D / "nhood_count.csv", index_col=0).values.astype(float)
z = (z + z.T) / 2
oe = np.log2(c / (np.outer(c.sum(1), c.sum(0)) / c.sum()))
n = len(z)

fig, axes = plt.subplots(1, 2, figsize=(15, 6.8), constrained_layout=True)
for ax, M, lim, title, lab, fmt in [
    (axes[0], z, 150, "Neighbourhood enrichment z-score", "z (colour clipped at ±150)", "{:.0f}"),
    (axes[1], oe, 3, "log2(observed / expected) neighbour pairs", "log2 O/E (clipped at ±3)", "{:.1f}"),
]:
    im = ax.imshow(np.clip(M, -lim, lim), cmap="RdBu_r", norm=TwoSlopeNorm(0, -lim, lim))
    for i in range(n):
        for j in range(n):
            v = M[i, j]
            ax.text(j, i, fmt.format(v), ha="center", va="center", fontsize=6.5,
                    color="white" if abs(v) > 0.6 * lim else "#222")
    for i in G:
        for j in G:
            ax.add_patch(Rectangle((j - .5, i - .5), 1, 1, fill=False, ec="black", lw=1.6))
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    lbl = [f"{k}*" if k in G else str(k) for k in range(n)]
    ax.set_xticklabels(lbl, fontsize=8); ax.set_yticklabels(lbl, fontsize=8)
    ax.set_xlabel("state"); ax.set_ylabel("state"); ax.set_title(title, fontsize=11)
    ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_visible(False)
    fig.colorbar(im, ax=ax, shrink=0.8, label=lab)
fig.suptitle("CRC P2 · Leiden start · K=16 · 6-NN graph · boxed = states 6, 8, 9, 10", fontsize=12)
for ext in ("png", "pdf"):
    fig.savefig(D / f"nhood_heatmap.{ext}", dpi=200)
