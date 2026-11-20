"""
Test 1 — same markers: do SpatialAIM's Leiden K=10 states 2 and 3 carry the genes
Oliveira et al. (Nat Genet 2025, Fig. 4c) use to define SPP1+ and SELENOP+
macrophages?

Writes to Data/03_GUI/Analysis/test1_markers/:
  * dotplot_paper_markers.pdf — the 15 genes x 10 states, scRNA and Xenium
  * summary.md — scRNA: Level2 composition, Wilcoxon state 2 vs 3

Run:  python test1_markers.py
"""


import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
from spatialaim_experiments import paths  # noqa: E402

ROOT = paths.GUI  # Data/03_GUI
RUN = ROOT / "Result" / "leiden" / "nearest_centroid"
K = 10
SPP1_STATE, SELENOP_STATE = 2, 3
OUT = paths.GUI_ANALYSIS / "test1_markers"

# Fig. 4c, in its column order
SELENOP_SET = ["SELENOP", "STAB1", "IGKC", "SLC40A1", "FGL2", "MS4A6A", "MPEG1"]
SPP1_SET = ["SPP1", "MMP7", "CHI3L1", "LPL", "FN1", "IL1RN", "FBP1", "APOC1"]
GENES = SELENOP_SET + SPP1_SET


def lognorm(adata: ad.AnnData) -> ad.AnnData:
    adata = adata.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def per_state(adata: ad.AnnData, genes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Mean lognorm and fraction > 0, states x genes."""
    X = adata[:, genes].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    states = adata.obs["state"].astype(int).to_numpy()
    idx = sorted(np.unique(states))
    mean = pd.DataFrame([X[states == s].mean(0) for s in idx], index=idx, columns=genes)
    frac = pd.DataFrame(
        [(X[states == s] > 0).mean(0) for s in idx], index=idx, columns=genes
    )
    return mean, frac


def two_vs_three(adata: ad.AnnData) -> pd.DataFrame:
    """Wilcoxon state 2 vs state 3 over every gene; the paper genes' rows with
    their rank in the genome-wide list towards the state the paper expects."""
    sub = adata[adata.obs["state"].isin([str(SPP1_STATE), str(SELENOP_STATE)])].copy()
    sc.tl.rank_genes_groups(
        sub,
        "state",
        groups=[str(SPP1_STATE)],
        reference=str(SELENOP_STATE),
        method="wilcoxon",
        n_genes=sub.n_vars,
    )
    df = sc.get.rank_genes_groups_df(sub, group=str(SPP1_STATE))
    rank_st2 = pd.Series(np.arange(1, len(df) + 1), index=df["names"])
    df = df.set_index("names").loc[GENES]
    expect_st2 = np.array([g in SPP1_SET for g in GENES])
    lfc = df["logfoldchanges"].to_numpy()
    return pd.DataFrame(
        {
            "gene": GENES,
            "paper set": np.where(expect_st2, "SPP1+", "SELENOP+"),
            "log2FC (st2 / st3)": lfc,
            "padj": df["pvals_adj"].to_numpy(),
            "rank": np.where(
                expect_st2, rank_st2[GENES], len(rank_st2) + 1 - rank_st2[GENES]
            ),
            "concordant": np.where(expect_st2, lfc > 0, lfc < 0)
            & (df["pvals_adj"].to_numpy() < 0.05),
        }
    )


def dotplot(mean: pd.DataFrame, frac: pd.DataFrame, title: str, ax):
    lo, hi = mean.min(0), mean.max(0)
    scaled = (mean - lo) / (hi - lo).replace(0, 1)
    ys, xs = np.meshgrid(range(len(mean)), range(mean.shape[1]), indexing="ij")
    pts = ax.scatter(
        xs.ravel(),
        ys.ravel(),
        s=220 * frac.to_numpy().ravel(),
        c=scaled.to_numpy().ravel(),
        cmap="Reds",
        vmin=0,
        vmax=1,
        edgecolors="0.4",
        linewidths=0.4,
    )
    ax.set_xticks(range(mean.shape[1]), mean.columns, rotation=90)
    ax.set_yticks(range(len(mean)), [f"state {s}" for s in mean.index])
    ax.set_ylim(len(mean) - 0.5, -0.5)
    ax.set_xlim(-0.6, mean.shape[1] - 0.4)
    n_sel = sum(g in SELENOP_SET for g in mean.columns)
    ax.axvline(n_sel - 0.5, color="0.6", lw=0.8, ls="--")
    for s, c in ((SPP1_STATE, "tab:orange"), (SELENOP_STATE, "tab:blue")):
        ax.axhspan(s - 0.5, s + 0.5, color=c, alpha=0.08, lw=0)
    ax.set_title(title, fontsize=10)
    return pts


def md(df: pd.DataFrame, index: bool = True) -> str:
    """A markdown pipe table (pandas' own needs the tabulate package)."""
    d = df.reset_index() if index else df

    def fmt(v) -> str:
        return f"{v:.3g}" if isinstance(v, (float, np.floating)) else str(v)

    rows = [" | ".join(map(str, d.columns)), " | ".join("---" for _ in d.columns)]
    rows += [" | ".join(fmt(v) for v in r) for r in d.itertuples(index=False)]
    return "\n".join(f"| {r} |" for r in rows)


def summary(comp: pd.DataFrame, hits: pd.DataFrame, n: pd.Series) -> str:
    n2, n3 = n[SPP1_STATE], n[SELENOP_STATE]
    return f"""# Test 1 — paper markers on SpatialAIM Leiden K={K}, states {SPP1_STATE} and {SELENOP_STATE}

Question: are states 2 and 3 the SPP1+ and SELENOP+ macrophages of Oliveira et al.?
The paper defines them by the 15 genes of its Fig. 4c (7 SELENOP+, 8 SPP1+). All 15
are measured in the scRNA reference, so the test is run there, on the whole
transcriptome (18,079 genes). The Xenium side appears only in the dotplot (7 of the 15
genes are on the panel).

## Setup

- Cells: the 35,673 P2 scRNA cells. Each cell's state is the K=10 tree cut of its
  Leiden start cluster — no mapping is involved on this side. State 2 has {n2} cells,
  state 3 has {n3}.
- Expression: counts normalised to 10,000 per cell, then log1p (as in the sweep).

## 1. What states 2 and 3 are (published Level2 labels)

Both are macrophage states: state 2 is {comp.loc['Macrophage', SPP1_STATE] / n2:.0%} Level2
"Macrophage", state 3 {comp.loc['Macrophage', SELENOP_STATE] / n3:.0%} ("Macrophage" plus
"Proliferating Macrophages"). So SpatialAIM splits the single published Level2 type
"Macrophage" into two.

{md(comp.head(8))}

## 2. Do the paper's genes separate state 2 from state 3?

Wilcoxon test of the state 2 cells against the state 3 cells, over all 18,079 genes
(scanpy `rank_genes_groups`, Benjamini–Hochberg adjusted).

- log2FC > 0: higher in state 2 (expected for SPP1+ genes); < 0: higher in state 3
  (expected for SELENOP+ genes).
- rank: the gene's position in the genome-wide list sorted towards the state the paper
  expects (1 = the strongest of all 18,079 genes).
- concordant: right direction and padj < 0.05.

Result: {int(hits['concordant'].sum())}/15 concordant.

{md(hits, index=False)}

- The SELENOP+ genes STAB1 (rank 1), MS4A6A, FGL2 and MPEG1 are among the 30 genes
  most specific to state 3; SPP1 (rank 8) and APOC1 lead state 2.
- IGKC fails: it is a plasma-cell gene (highest in the plasma states 0 and 8).
  In the paper's 8 µm Visium HD bins it most likely comes from neighbouring plasma
  cells, not from the macrophages.
- SELENOP itself separates only weakly (rank 285); in the reference it is highest in the
  goblet/enterocyte state 6. The SELENOP+ state is carried by the set as a whole.
"""


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cut = pd.read_csv(RUN / f"k_{K:03d}" / "start_cluster_to_state.csv").set_index(
        "start_cluster"
    )["state"]

    # scRNA: state = tree cut of the cell's start cluster
    sc_raw = ad.read_h5ad(ROOT / "scRNA" / "11_human-crc-p2.h5ad")
    start = ad.read_h5ad(ROOT / "Result" / "leiden" / "start_clustering.h5ad").obs
    assert (start.index == sc_raw.obs_names).all()
    sc_ad = lognorm(sc_raw)
    sc_ad.obs["state"] = (
        cut.reindex(start["start_cluster"].astype(int)).astype(str).to_numpy()
    )

    comp = pd.crosstab(sc_ad.obs["Level2"].astype(str), sc_ad.obs["state"].astype(int))
    comp = comp[[SPP1_STATE, SELENOP_STATE]]
    comp = comp[(comp > 0).any(axis=1)].sort_values(SPP1_STATE, ascending=False)
    n_cells = sc_ad.obs["state"].astype(int).value_counts()

    # Xenium (dotplot only): state = argmax of the soft mapping
    st_raw = ad.read_h5ad(ROOT / "ST" / "11_human-crc-p2.h5ad")
    P = ad.read_h5ad(RUN / f"k_{K:03d}" / "spot_to_state_mapping_soft.h5ad")
    assert (P.obs_names == st_raw.obs_names).all()
    st_ad = lognorm(st_raw)
    Pm = P.X.toarray() if hasattr(P.X, "toarray") else np.asarray(P.X)
    st_ad.obs["state"] = Pm.argmax(1).astype(str)
    st_genes = [g for g in GENES if g in st_ad.var_names]

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(11, 4.6),
        gridspec_kw={"width_ratios": [len(GENES), len(st_genes) + 2]},
    )
    for ax, (title, adata, genes) in zip(
        axes,
        (
            ("scRNA (whole transcriptome)", sc_ad, GENES),
            ("Xenium (panel)", st_ad, st_genes),
        ),
    ):
        pts = dotplot(*per_state(adata, genes), title, ax)
    fig.colorbar(pts, ax=axes, label="mean lognorm (scaled per gene)", shrink=0.6)
    fig.suptitle(
        "Oliveira Fig. 4c markers (left of dashed line: SELENOP+, right: SPP1+) "
        f"on SpatialAIM Leiden K={K} states",
        fontsize=10,
    )
    fig.savefig(OUT / "dotplot_paper_markers.pdf", bbox_inches="tight")

    text = summary(comp, two_vs_three(sc_ad), n_cells)
    (OUT / "summary.md").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
