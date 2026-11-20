"""
Test 2 — same niches: do SpatialAIM's Xenium cells of Leiden K=10 state 2 / 3 sit
where Oliveira et al. labelled Visium HD bins SPP1+ / SELENOP+ macrophages?

Visium HD (authors' P2CRC_Metadata.parquet, GitHub MetaData) and Xenium are adjacent
5 um sections in their own frames, so:
  1. Register: rigid (rotation, flip, translation) fit of the Visium HD bins onto the
     Xenium cells, maximising the overlap of tumour (Visium "Tumor" bins / our tumour
     state 4) and non-tumour tissue maps. No macrophage label is used.
  2. Restrict Xenium state 2/3 cells to the paper's region: the nearest registered
     bin (within 8 um) is one of their "50 micron" tumour-periphery bins.
  3. Per SPP1+ / SELENOP+ bin: our state 2 and 3 cells (of step 2) within 24 um of its
     centre; which state is in the majority, per bin label.

Writes to Data/03_GUI/Analysis/test2_markers/: registration.pdf, densities.pdf (counts per 100 um tile),
summary.md, transform.json.
Run:  python test2_markers.py
"""

import json

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import ndimage, signal, sparse  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402
from spatialaim_experiments import paths  # noqa: E402

ROOT = paths.GUI  # Data/03_GUI
RUN = ROOT / "Result" / "leiden" / "nearest_centroid"
META = ROOT / "_downloads" / "visiumhd_p2_metadata" / "P2CRC_Metadata.parquet"
OUT = paths.GUI_ANALYSIS / "test2_markers"
K = 10
SPP1_STATE, SELENOP_STATE, TUMOUR_STATE = 2, 3, 4
BIN_UM = 8.0
PX = 40.0  # registration raster (um)
CANVAS = 12000.0  # registration canvas side (um), centred on each dataset
R_UM = 24.0  # cells counted around each bin centre
TILE_UM = 100.0  # density heatmap tiles
N_PERM = 1000  # state 2/3 label shuffles
RNG = np.random.default_rng(0)


# --------------------------------------------------------------------------- #
# Data
# --------------------------------------------------------------------------- #
def load_visium() -> pd.DataFrame:
    m = pd.read_parquet(
        META, columns=["barcode", "Periphery", "MacrophageSubtype"]
    )
    rc = m["barcode"].str.extract(r"s_008um_(\d+)_(\d+)-1").astype(int)
    # Bin grid -> um (row = y, col = x); the X/Y columns are in arbitrary units.
    m["x"], m["y"] = rc[1] * BIN_UM, rc[0] * BIN_UM
    return m[m["Periphery"].notna()].reset_index(drop=True)  # in-tissue bins


def load_xenium() -> tuple[np.ndarray, np.ndarray]:
    st = ad.read_h5ad(ROOT / "ST" / "11_human-crc-p2.h5ad", backed="r")
    xy = np.asarray(st.obsm["spatial"][:], dtype=float)
    P = ad.read_h5ad(RUN / f"k_{K:03d}" / "spot_to_state_mapping_soft.h5ad")
    assert (P.obs_names == st.obs_names).all()
    Pm = P.X.toarray() if hasattr(P.X, "toarray") else np.asarray(P.X)
    return xy, Pm.argmax(1)


# --------------------------------------------------------------------------- #
# Registration
# --------------------------------------------------------------------------- #
def raster(xy: np.ndarray, centre: np.ndarray) -> np.ndarray:
    n = int(CANVAS / PX)
    edges = np.linspace(-CANVAS / 2, CANVAS / 2, n + 1)
    img, *_ = np.histogram2d(xy[:, 1] - centre[1], xy[:, 0] - centre[0], [edges, edges])
    img = ndimage.gaussian_filter(img, 1.0)
    return img / (img.max() or 1.0)


def channels(xy: np.ndarray, tumour: np.ndarray, centre: np.ndarray) -> list[np.ndarray]:
    imgs = [raster(xy[tumour], centre), raster(xy[~tumour], centre)]
    return [i - i.mean() for i in imgs]


def transform(xy, centre_from, angle, flip, centre_to, shift):
    p = xy - centre_from
    if flip:
        p = p * np.array([-1.0, 1.0])
    c, s = np.cos(angle), np.sin(angle)
    return p @ np.array([[c, s], [-s, c]]) + centre_to + shift


def best_shift(moving: list[np.ndarray], fixed: list[np.ndarray]) -> tuple[float, np.ndarray]:
    """Translation (um) maximising the summed cross-correlation of the channels."""
    cc = sum(
        signal.fftconvolve(f, m[::-1, ::-1], mode="same") for m, f in zip(moving, fixed)
    )
    iy, ix = np.unravel_index(np.argmax(cc), cc.shape)
    mid = np.array(cc.shape) // 2
    norm = np.sqrt(sum((m**2).sum() for m in moving) * sum((f**2).sum() for f in fixed))
    return cc[iy, ix] / norm, np.array([ix - mid[1], iy - mid[0]]) * PX


def register(vis: pd.DataFrame, xen_xy: np.ndarray, xen_tumour: np.ndarray):
    vxy = vis[["x", "y"]].to_numpy()
    vt = (vis["Periphery"] == "Tumor").to_numpy()
    cv, cx = vxy.mean(0), xen_xy.mean(0)
    fixed = channels(xen_xy, xen_tumour, cx)
    best = (-np.inf, None)
    for flip in (False, True):
        for deg in np.arange(0, 360, 2.0):
            moved = transform(vxy, cv, np.deg2rad(deg), flip, cx, 0.0)
            score, shift = best_shift(channels(moved, vt, cx), fixed)
            if score > best[0]:
                best = (score, (deg, flip, shift))
    # Refine the angle around the coarse optimum.
    score, (deg0, flip, shift) = best
    for deg in np.arange(deg0 - 2, deg0 + 2.01, 0.25):
        moved = transform(vxy, cv, np.deg2rad(deg), flip, cx, 0.0)
        s, sh = best_shift(channels(moved, vt, cx), fixed)
        if s > score:
            score, deg0, shift = s, deg, sh
    params = dict(
        angle_deg=float(deg0),
        flip=bool(flip),
        shift_um=shift,
        score=float(score),
        centre_visium_um=cv,
        centre_xenium_um=cx,
    )
    reg = transform(vxy, cv, np.deg2rad(deg0), flip, cx, shift)
    return reg, params


def tumour_dice(reg: np.ndarray, vt: np.ndarray, xen_xy, xen_tumour) -> float:
    """Dice of the tumour masks (tumour share > 0.5 per 100 um pixel) where both
    sections have tissue."""
    lo = np.minimum(reg.min(0), xen_xy.min(0))
    hi = np.maximum(reg.max(0), xen_xy.max(0))
    edges = [np.arange(lo[i], hi[i] + 100, 100) for i in (0, 1)]

    def share(xy, t):
        a, *_ = np.histogram2d(xy[:, 0], xy[:, 1], edges)
        b, *_ = np.histogram2d(xy[t, 0], xy[t, 1], edges)
        return a > 0, np.where(a > 0, b / np.maximum(a, 1), 0) > 0.5

    tv, mv = share(reg, vt)
    tx, mx = share(xen_xy, xen_tumour)
    both = tv & tx
    a, b = mv & both, mx & both
    return 2 * (a & b).sum() / (a.sum() + b.sum())


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    vis = load_visium()
    xen_xy, hard = load_xenium()
    xen_tumour = hard == TUMOUR_STATE

    reg, params = register(vis, xen_xy, xen_tumour)
    vt = (vis["Periphery"] == "Tumor").to_numpy()
    dice = tumour_dice(reg, vt, xen_xy, xen_tumour)
    write_transform(params, dice)

    # Xenium state 2/3 cells in the paper's region: the nearest registered in-tissue
    # bin (within one bin width) is one of their "50 micron" periphery bins.
    mac = np.isin(hard, [SPP1_STATE, SELENOP_STATE])
    dist, idx = cKDTree(reg).query(xen_xy[mac], distance_upper_bound=BIN_UM)
    found = np.isfinite(dist)
    periphery = (vis["Periphery"] == "50 micron").to_numpy()
    region = np.zeros(len(idx), dtype=bool)
    region[found] = periphery[idx[found]]
    cells_xy, state = xen_xy[mac][region], hard[mac][region]

    sub = vis["MacrophageSubtype"].notna().to_numpy()
    bins_xy = reg[sub]
    bins_spp1 = (vis.loc[sub, "MacrophageSubtype"] == "Macrophage-SPP1+").to_numpy()

    near = near_matrix(cells_xy, bins_xy)
    per_bin = bin_votes(near, state, bins_spp1)
    perm = permutation_test(near, state, bins_spp1)

    plot_registration(reg, vt, xen_xy, xen_tumour, bins_xy, bins_spp1, cells_xy, state)
    plot_densities(xen_xy, bins_xy, bins_spp1, cells_xy, state)
    write_summary(params, dice, vis, sub, bins_spp1, mac, region, state, per_bin, perm)


def near_matrix(cells_xy, bins_xy) -> sparse.csr_matrix:
    """bins x cells, 1 where the cell's centre lies within R_UM of the bin's centre."""
    hits = cKDTree(cells_xy).query_ball_point(bins_xy, R_UM)
    rows = np.repeat(np.arange(len(hits)), [len(h) for h in hits])
    cols = np.concatenate([np.asarray(h, dtype=int) for h in hits])
    return sparse.csr_matrix(
        (np.ones(len(cols)), (rows, cols)), shape=(len(bins_xy), len(cells_xy))
    )


def _shares(near, state, bins_spp1) -> dict[str, np.ndarray]:
    """Per bin label (SPP1+, SELENOP+): shares of its bins with no state 2/3 cell,
    more state 2, more state 3, a tie."""
    n2 = near @ (state == SPP1_STATE)
    n3 = near @ (state == SELENOP_STATE)
    out = {k: np.zeros(2) for k in ("none", "more2", "more3", "tie")}
    for i, mask in enumerate((bins_spp1, ~bins_spp1)):
        a, b = n2[mask], n3[mask]
        out["none"][i] = ((a + b) == 0).mean()
        out["more2"][i] = (a > b).mean()
        out["more3"][i] = (b > a).mean()
        out["tie"][i] = ((a == b) & (a > 0)).mean()
    return out


def bin_votes(near, state, bins_spp1) -> pd.DataFrame:
    s = _shares(near, state, bins_spp1)
    return pd.DataFrame(
        {
            "bin label": ["SPP1+", "SELENOP+"],
            "bins": [int(bins_spp1.sum()), int((~bins_spp1).sum())],
            "no state 2/3 cell": s["none"],
            "more state 2": s["more2"],
            "more state 3": s["more3"],
            "tie": s["tie"],
        }
    )


def permutation_test(near, state, bins_spp1) -> pd.DataFrame:
    """Shuffle the state 2/3 labels among the periphery cells (positions and the
    number of each state fixed) N_PERM times. Tested per bin label: the matching
    share, "more state 2" for SPP1+ bins and "more state 3" for SELENOP+ bins."""
    obs = _shares(near, state, bins_spp1)
    null = [_shares(near, RNG.permutation(state), bins_spp1) for _ in range(N_PERM)]
    rows = []
    for i, (label, key, col) in enumerate(
        (("SPP1+", "more2", "more state 2"), ("SELENOP+", "more3", "more state 3"))
    ):
        v = np.array([n[key][i] for n in null])
        rows.append(
            {
                "bin label": label,
                "tested share": col,
                "observed": obs[key][i],
                "p": (1 + (v >= obs[key][i]).sum()) / (1 + N_PERM),
            }
        )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #
def plot_registration(reg, vt, xen_xy, xen_tumour, bins_xy, bins_spp1, cells_xy, state):
    fig, axes = plt.subplots(1, 2, figsize=(13, 6.2))
    ax = axes[0]
    take = RNG.choice(len(xen_xy), min(len(xen_xy), 150_000), replace=False)
    ax.scatter(*xen_xy[take].T, s=0.2, c=np.where(xen_tumour[take], "0.35", "0.85"),
               rasterized=True, lw=0)
    ax.scatter(*reg[vt][:: 3].T, s=0.2, c="tab:red", alpha=0.15, rasterized=True, lw=0)
    ax.set_title("Registration: Xenium cells (dark = tumour state 4)\n"
                 "+ Visium HD tumour bins (red)", fontsize=9)
    ax = axes[1]
    ax.scatter(*xen_xy[take].T, s=0.2, c="0.9", rasterized=True, lw=0)
    for s, c in ((SELENOP_STATE, "tab:blue"), (SPP1_STATE, "tab:orange")):
        ax.scatter(*cells_xy[state == s].T, s=1.5, c=c, lw=0, rasterized=True,
                   label=f"Xenium state {s}")
    ax.scatter(*bins_xy[~bins_spp1].T, s=9, marker="s", facecolors="none",
               edgecolors="navy", lw=0.6, label="Visium HD SELENOP+ bin")
    ax.scatter(*bins_xy[bins_spp1].T, s=9, marker="s", facecolors="none",
               edgecolors="darkred", lw=0.6, label="Visium HD SPP1+ bin")
    ax.legend(markerscale=3, fontsize=7, loc="lower left")
    ax.set_title("Their 50 um tumour periphery: our states 2/3 vs their subtype bins",
                 fontsize=9)
    for a in axes:
        a.set_aspect("equal")
        a.invert_yaxis()
        a.set_xlabel("x (um)")
    axes[0].set_ylabel("y (um)")
    fig.savefig(OUT / "registration.pdf", bbox_inches="tight", dpi=200)


def plot_densities(xen_xy, bins_xy, bins_spp1, cells_xy, state) -> None:
    """Counts per TILE_UM x TILE_UM tile on one grid over the Xenium section: their
    SPP1+ / SELENOP+ bins (registered) and our periphery state 2 / 3 cells."""
    lo = np.floor(xen_xy.min(0) / TILE_UM) * TILE_UM
    hi = np.ceil(xen_xy.max(0) / TILE_UM) * TILE_UM
    ex, ey = (np.arange(lo[i], hi[i] + TILE_UM, TILE_UM) for i in (0, 1))
    panels = (
        ("Visium HD SPP1+ bins", bins_xy[bins_spp1], "Reds"),
        ("Visium HD SELENOP+ bins", bins_xy[~bins_spp1], "Blues"),
        (f"Xenium state {SPP1_STATE} cells", cells_xy[state == SPP1_STATE], "Reds"),
        (f"Xenium state {SELENOP_STATE} cells", cells_xy[state == SELENOP_STATE], "Blues"),
    )
    fig, axes = plt.subplots(2, 2, figsize=(12, 12.5))
    for ax, (title, xy, cmap) in zip(axes.ravel(), panels):
        counts, *_ = np.histogram2d(xy[:, 0], xy[:, 1], [ex, ey])
        im = ax.pcolormesh(
            ex, ey, np.ma.masked_equal(counts.T, 0), cmap=cmap, vmin=0, shading="flat"
        )
        fig.colorbar(im, ax=ax, shrink=0.75, label=f"count per {TILE_UM:.0f} um tile")
        ax.set_title(f"{title} (n = {len(xy):,})", fontsize=10)
        ax.set_xlim(ex[0], ex[-1])
        ax.set_ylim(ey[-1], ey[0])  # y down, as in registration.pdf
        ax.set_aspect("equal")
        ax.set_facecolor("0.96")
        ax.set_xlabel("x (um)")
        ax.set_ylabel("y (um)")
    fig.suptitle(
        f"Counts per {TILE_UM:.0f} x {TILE_UM:.0f} um tile, Xenium frame after "
        "registration; cells: the periphery set of section 2",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(OUT / "densities.pdf", bbox_inches="tight")


def write_transform(params: dict, dice: float) -> None:
    """The rigid Visium HD -> Xenium transform, enough to reproduce it."""
    a = np.deg2rad(params["angle_deg"])
    doc = {
        "direction": "Visium HD bins (moving) -> Xenium cells (fixed); Xenium untouched",
        "visium_input": "bin grid in um: x = col * 8, y = row * 8, "
        "from barcode s_008um_<row>_<col>-1",
        "xenium_frame": "obsm['spatial'] of ST/11_human-crc-p2.h5ad (x_centroid, "
        "y_centroid, um)",
        "formula": "p = [x, y] - centre_visium_um; if flip: p = p * [-1, 1]; "
        "out = p @ [[cos a, sin a], [-sin a, cos a]] + centre_xenium_um + shift_um",
        "angle_deg": params["angle_deg"],
        "flip": params["flip"],
        "centre_visium_um": [float(v) for v in params["centre_visium_um"]],
        "centre_xenium_um": [float(v) for v in params["centre_xenium_um"]],
        "shift_um": [float(v) for v in params["shift_um"]],
        "rotation_matrix": [[np.cos(a), np.sin(a)], [-np.sin(a), np.cos(a)]],
        "fit": {
            "normalised_cross_correlation": params["score"],
            "tumour_mask_dice": float(dice),
            "raster_um": PX,
        },
    }
    (OUT / "transform.json").write_text(json.dumps(doc, indent=2), encoding="utf-8")


def md(df: pd.DataFrame, index: bool = True) -> str:
    d = df.reset_index() if index else df

    def fmt(v) -> str:
        return f"{v:.3g}" if isinstance(v, (float, np.floating)) else str(v)

    rows = [" | ".join(map(str, d.columns)), " | ".join("---" for _ in d.columns)]
    rows += [" | ".join(fmt(v) for v in r) for r in d.itertuples(index=False)]
    return "\n".join(f"| {r} |" for r in rows)


def write_summary(params, dice, vis, sub, bins_spp1, mac, region, state, per_bin, perm):
    per_bin, perm = per_bin.copy(), perm.copy()
    for c in per_bin.columns[2:]:
        per_bin[c] = [f"{v:.0%}" for v in per_bin[c]]
    for c in perm.columns[2:-1]:
        perm[c] = [f"{v:.0%}" for v in perm[c]]
    perm["p"] = [f"{v:.3g}" for v in perm["p"]]
    base = (state == SPP1_STATE).mean()
    sh = params["shift_um"]
    text = f"""# Test 2 — our Xenium states 2/3 vs their Visium HD SPP1+/SELENOP+ bins

Data: Visium HD P2CRC 8 um bins with the authors' labels (`P2CRC_Metadata.parquet`,
GitHub `10XGenomics/HumanColonCancer_VisiumHD/MetaData`): {int(sub.sum())} macrophage
bins carry a subtype, {int(bins_spp1.sum())} SPP1+ and {int((~bins_spp1).sum())} SELENOP+, all in
their 50 um tumour periphery. Xenium: SpatialAIM Leiden K={K}, hard assignment.
The two assays ran on adjacent 5 um sections (Methods), so the same niches, not the
same cells.

## 1. Registration

Rigid fit (rotation, mirror, translation; no scaling) of the {len(vis):,} in-tissue bins
onto the Xenium cells, maximising the cross-correlation of two maps on a {PX:.0f} um grid:
tumour (Visium "Tumor" bins / Xenium state {TUMOUR_STATE}, 92 % Tumor in the reference) and
the remaining tissue. No macrophage label enters the fit.

- rotation {params['angle_deg']:.2f}°, mirrored: {params['flip']}, translation
  ({sh[0]:.0f}, {sh[1]:.0f}) um; normalised correlation {params['score']:.2f}
- tumour-mask Dice (100 um pixels, tumour share > 0.5, where both have tissue):
  {dice:.2f}

See `registration.pdf` (left).

## 2. Cells compared

Xenium state {SPP1_STATE}/{SELENOP_STATE} cells: {int(mac.sum()):,}; of these {int(region.sum()):,} lie in the
paper's region, i.e. their nearest registered bin (within {BIN_UM:.0f} um) is one of the
authors' "50 micron" tumour-periphery bins (the ring outside their "Tumor" bins, so
macrophages inside the tumour are excluded): {int((state == SPP1_STATE).sum()):,} in state {SPP1_STATE},
{int((state == SELENOP_STATE).sum()):,} in state {SELENOP_STATE}.

## 3. Per bin, our cells

Per SPP1+ or SELENOP+ bin, count our state {SPP1_STATE} and state {SELENOP_STATE} cells (of the
{len(state):,} periphery cells above) whose centres lie within {R_UM:.0f} um of the bin's centre.
Per bin label, the shares of its bins (summing to 100 %):

- no state 2/3 cell: no state 2 or 3 cell within {R_UM:.0f} um
- more state 2 / more state 3: more cells of that state
- tie: as many state 2 as state 3 cells (at least one each)

State {SPP1_STATE} makes up only {base:.0%} of the periphery state 2/3 cells, so a bin can be
"more state 3" by abundance alone; section 4 accounts for that.

{md(per_bin, index=False)}

## 4. Permutation test

Shuffle the state {SPP1_STATE}/{SELENOP_STATE} labels among the {len(state):,} periphery cells ({N_PERM:,} times):
positions and the number of cells per state stay fixed, only which cell carries which
state changes. Recompute section 3 each time. Tested per bin label: the matching share
("more state {SPP1_STATE}" for SPP1+ bins, "more state {SELENOP_STATE}" for SELENOP+ bins);
p = share of shuffles reaching the observed value (one-sided, smallest possible
{1 / (1 + N_PERM):.3g}).

{md(perm, index=False)}

The shuffle treats cells as independent; same-state cells cluster in space, so the
p-values are optimistic.
"""
    (OUT / "summary.md").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
