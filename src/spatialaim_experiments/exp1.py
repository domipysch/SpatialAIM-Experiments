"""Experiment 1's results table, shared by the reconstruction and biology figures
(sections 03 and 04).

``summary/best_per_run.csv`` is built on the cluster by ``experiment1/04_summarize.py`` and
pulled to ``Data/02_Results/experiment1/summary/`` (``paths.EXP1_SUMMARY``). Each
row is one stage-B run at the one K the figures compare:

    baseline     annotation start, K = T (the identity cut: every spot maps
                 straight onto the author's annotation)
    spatialaim   Leiden start, K = best.overall from k_selection.json

Both are read from the same column of the same metric table, so the comparison
is apples to apples.
"""

from __future__ import annotations
import pandas as pd
from . import paths


def load() -> pd.DataFrame:
    """``summary/best_per_run.csv``, one row per stage-B run (see module docstring)."""
    path = paths.EXP1_SUMMARY / "best_per_run.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} missing - pull the Experiment 1 summary (experiment1/README.md)")
    return pd.read_csv(path, keep_default_na=False, na_values=[""])


def resources(stage: str) -> pd.DataFrame:
    """``summary/resources_stage_<stage>.csv`` (stage "a" = scaffold builds, "b" =
    mapping runs): time and memory per run, from its resources.json."""
    return pd.read_csv(paths.EXP1_SUMMARY / f"resources_stage_{stage}.csv",
                       keep_default_na=False, na_values=[""])


def paired(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    """Put each SpatialAIM run next to the annotation run it is compared with.

    Input
        df      best_per_run.csv: one row per run, role "spatialaim"
                (Leiden start, proposed K) or "baseline" (annotation start, K = T)
        metric  the metric to compare, e.g. "cossim_hard_raw_spot"

    Output
        one row per (pair, tree, mapper), i.e. per SpatialAIM run, with
          all columns of the SpatialAIM row (``metric``, ``k``, ...)
          <metric>_base, k_base, ann_key_base   the same for the baseline
          delta                                 metric - metric_base

    Both sides use the same pair, tree and mapper, so the difference is only the
    starting partition. Where a reference has two annotations, the baseline is
    the one with the higher ``score_overall`` -- the score SpatialAIM picks its
    own K by. Choosing it on the plotted metric would favour the baseline per
    panel; always taking the first annotation would be arbitrary, and it matters
    (on 01_mouse-ssp SpatialAIM beats ``cellType``, 3 types, by +0.26 per-spot
    cosine and loses to ``cellTypeMinor``, 22 types, by -0.05).
    """
    join = ["pair_dir", "tree", "mapper"]
    spatialaim = df[df["role"] == "spatialaim"]
    base = df[df["role"] == "baseline"]
    base = base.loc[base.groupby(join)["score_overall"].idxmax()]
    out = spatialaim.merge(base[join + [metric, "k", "ann_key"]],
                           on=join, suffixes=("", "_base"), validate="one_to_one")
    out["delta"] = out[metric] - out[f"{metric}_base"]
    return out
