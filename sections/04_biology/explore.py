"""Results 2.4 -- biological and topological plausibility of AIM's cell states.

Numbers only; the figure lives in ``panels.py``. Reuses the Results 2.3 table
(``reconstruction/load.py``), which now carries all four non-reconstruction
metrics per run root.

The comparison this section makes, fixed by the section plan:

    AIM        Leiden over-clustering, **average** linkage, K = best.overall
    annotation the author's annotation mapped as it stands (K = T), and where a
               reference carries two annotations, **the one whose identity cut
               scores highest overall** -- the same harmonic mean across the
               three criteria by which AIM picks its own K. Choosing it on the
               measure a given panel plots would hand the baseline a per-panel
               advantage AIM is never given.
    random     0 for all four metrics, by construction: the two spatial metrics
               are z-scores against a spot-label shuffle, and modularity is
               defined as the excess over a degree-preserving random partition.

Both sides use the same mapper, so what differs is only the starting partition.
"""

from __future__ import annotations

import pandas as pd

from aim_experiments import exp1 as L

METRICS = {
    "nhood_mean_self_zscore": "neighbourhood enrichment z",
    "local_purity_zscore": "local purity z",
    "modularity_st_expression": "ST-expression modularity",
    "modularity_shared": "SC shared-gene modularity",
}
SPATIAL = ["nhood_mean_self_zscore", "local_purity_zscore"]
COHERENCE = ["modularity_st_expression", "modularity_shared"]

TREE = "average"


def frames(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    """One row per (pair, mapper): AIM against its annotation baseline."""
    aim = df[(df["role"] == "aim") & (df["tree"] == TREE)]
    base = df[(df["role"] == "baseline") & (df["tree"] == TREE)]
    join = ["pair_dir", "mapper"]
    base = base.loc[base.groupby(join)["score_overall"].idxmax()]
    out = aim.merge(
        base[join + [metric, "k", "ann_key"]],
        on=join, suffixes=("", "_base"), validate="one_to_one",
    )
    out["delta"] = out[metric] - out[f"{metric}_base"]
    out["dk"] = out["k"] - out["k_base"]
    return out


def check_tree_independence(df: pd.DataFrame) -> None:
    """At K = T the tree cut is the identity, so ward and average must agree."""
    base = df[df["role"] == "baseline"]
    wide = base.pivot_table(
        index=["pair_dir", "mapper", "ann_key"], columns="tree",
        values=list(METRICS) + list(L.METRICS),
    )
    worst = 0.0
    for metric in list(METRICS) + list(L.METRICS):
        d = (wide[(metric, "ward")] - wide[(metric, "average")]).abs()
        worst = max(worst, float(d.max()))
    print(f"baseline ward vs average, largest absolute difference over all "
          f"six metrics: {worst:.3g}")


def main() -> None:
    df = L.load()
    print(f"{len(df)} run roots, {df['pair_dir'].nunique()} pairs, "
          f"{df['reference'].nunique()} references\n")
    check_tree_independence(df)

    aim = df[(df["role"] == "aim") & (df["tree"] == TREE)]
    print(f"\nAIM side: {len(aim)} runs "
          f"({aim['pair_dir'].nunique()} pairs x {aim['mapper'].nunique()} mappers)")
    print("AIM K, per reference:")
    print(aim.groupby("reference")["k"].agg(["min", "median", "max"]).to_string())

    for metric, label in METRICS.items():
        pr = frames(df, metric)
        print(f"\n{'=' * 72}\n{label}   [{metric}]\n{'=' * 72}")
        print(f"{len(pr)} comparisons")

        # --- vs random -------------------------------------------------------
        print("\n  vs random (chance = 0)")
        print(f"    AIM         mean {pr[metric].mean():9.4f}   "
              f"above 0: {(pr[metric] > 0).mean():6.1%}")
        print(f"    annotation  mean {pr[metric + '_base'].mean():9.4f}   "
              f"above 0: {(pr[metric + '_base'] > 0).mean():6.1%}")

        # --- vs the annotation ----------------------------------------------
        per_ref = pr.groupby("reference").agg(
            aim=(metric, "mean"), ann=(f"{metric}_base", "mean"),
            delta=("delta", "mean"), n=("delta", "size"),
            k=("k", "mean"), t=("k_base", "mean"),
        )
        print("\n  vs annotation (AIM - annotation), per reference")
        print(per_ref.round(4).to_string())
        print(f"    mean over references {per_ref['delta'].mean():+.4f}   "
              f"references positive {int((per_ref['delta'] > 0).sum())}/{len(per_ref)}   "
              f"comparisons positive {(pr['delta'] > 0).mean():.1%}")

        print("\n  per mapper (mean over references, so retina cannot dominate)")
        by = (pr.groupby(["reference", "mapper"])
                .agg(aim=(metric, "mean"), ann=(f"{metric}_base", "mean"),
                     delta=("delta", "mean"))
                .groupby("mapper").mean())
        print(by.round(4).to_string())

    # --- how the two criteria move together ---------------------------------
    print(f"\n{'=' * 72}\ncorrelation of the four deltas (per reference means)\n{'=' * 72}")
    cols = {}
    for metric in METRICS:
        cols[metric] = frames(df, metric).groupby(["reference", "mapper"])["delta"].mean()
    print(pd.DataFrame(cols).corr(method="spearman").round(3).to_string())


if __name__ == "__main__":
    main()
