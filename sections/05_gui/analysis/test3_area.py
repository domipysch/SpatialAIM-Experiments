"""
Test 3 — cell size: do SpatialAIM's Leiden K=10 states 2 (SPP1+) and 3 (SELENOP+)
differ in segmented Xenium cell area?

Writes to Data/03_GUI/Analysis/test3_area/:
  * summary.md — per-state areas, median difference, Mann-Whitney U

Run:  python test3_area.py
"""


import anndata as ad
import numpy as np
import pandas as pd
import tifffile
from scipy import stats
from skimage.measure import regionprops_table
from spatialaim_experiments import paths

ROOT = paths.GUI  # Data/03_GUI
RUN = ROOT / "Result" / "leiden" / "nearest_centroid"
ST = ROOT / "ST" / "11_human-crc-p2.h5ad"
MASK = ROOT / "ST" / "Segmentation" / "11_human-crc-p2.tif"
K = 10
SPP1_STATE, SELENOP_STATE = "state_2", "state_3"
OUT = paths.GUI_ANALYSIS / "test3_area"
RNG = np.random.default_rng(0)


def load() -> tuple[pd.DataFrame, float]:
    """Hard state and area (um^2) per Xenium cell; the mask's pixel size."""
    st = ad.read_h5ad(ST, backed="r")
    seg = st.obs["segmentation_id"].to_numpy().astype(np.int64)
    pixel_um = float(st.uns["segmentation_pixel_size_um"])
    m = ad.read_h5ad(RUN / f"k_{K:03d}" / "spot_to_state_mapping_soft.h5ad")
    assert (m.obs_names == st.obs_names).all()
    state = np.array(m.var_names)[np.asarray(m.X).argmax(1)]
    mask = tifffile.imread(MASK)
    props = pd.DataFrame(regionprops_table(mask, properties=("label", "area")))
    area = props.set_index("label")["area"].reindex(seg).to_numpy(float) * pixel_um**2
    area[seg == 0] = np.nan
    return pd.DataFrame({"state": state, "area": area}).dropna(), pixel_um


def main() -> None:
    OUT.mkdir(exist_ok=True)
    df, pixel_um = load()
    a2 = df.area[df.state == SPP1_STATE].to_numpy()
    a3 = df.area[df.state == SELENOP_STATE].to_numpy()

    boot = [
        np.median(RNG.choice(a2, len(a2))) - np.median(RNG.choice(a3, len(a3)))
        for _ in range(2000)
    ]
    ci = np.percentile(boot, [2.5, 97.5])
    u = stats.mannwhitneyu(a2, a3, alternative="two-sided")
    auc = u.statistic / (len(a2) * len(a3))

    rows = "\n".join(
        f"| {s} | {len(a)} | {np.percentile(a, 25):.1f} | {np.median(a):.1f} | "
        f"{np.percentile(a, 75):.1f} | {a.mean():.1f} |"
        for s, a in ((SPP1_STATE, a2), (SELENOP_STATE, a3))
    )

    md = f"""# Test 3 — cell area of SpatialAIM Leiden K={K} states 2 (SPP1+) and 3 (SELENOP+)

Question: do the two macrophage states differ in size in the Xenium data?

## Setup

- Cells: all {len(df):,} Xenium cells, hard state = argmax of the K={K} nearest-centroid mapping.
- Area: segmentation-mask pixels per cell (skimage `regionprops_table`, as SpatialAIM's shape
  features) x {pixel_um}^2 um^2.

## 1. Area per state (um^2)

| state | cells | Q1 | median | Q3 | mean |
| --- | --- | --- | --- | --- | --- |
{rows}

## 2. State 2 vs state 3

- Median {np.median(a2):.1f} vs {np.median(a3):.1f} um^2 (ratio {np.median(a2) / np.median(a3):.2f}); median difference
  {np.median(a2) - np.median(a3):.1f} um^2, bootstrap 95 % CI [{ci[0]:.1f}, {ci[1]:.1f}].
- Mann-Whitney U (two-sided): p = {u.pvalue:.2g}; P(state 2 > state 3) = {auc:.2f}.


## Caveats

- Segmented area, not true cell size; segmentation quality can differ between niches.
"""
    (OUT / "summary.md").write_text(md, encoding="utf-8")
    print(f"wrote {OUT / 'summary.md'}")


if __name__ == "__main__":
    main()
