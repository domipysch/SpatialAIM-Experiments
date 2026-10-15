"""Experiment 1's results table, shared by the reconstruction and biology figures
(sections 02 and 04).

``summary/best_per_run.csv`` is built on the cluster by ``experiment1/04_summarize.py`` and
pulled to ``Data/02_Results/experiment1/summary/`` (``paths.EXP1_SUMMARY``). Each
row is one stage-B run at the one K the figures compare:

    baseline     annotation start, K = T (the identity cut: every spot maps
                 straight onto the author's annotation)
    spatialaim        Leiden start, K = best.overall from k_selection.json
    annotation_best   annotation start, K = best.overall

The figures compare against the baseline the better of the two SpatialAIM starts
(``best_start``). All are read from the same column of the same metric table, so
the comparison is apples to apples.
"""

from __future__ import annotations
import logging
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


def annotation_best(df: pd.DataFrame) -> pd.DataFrame:
    """The annotation-start runs at K = best.overall: the ``annotation_best`` rows
    of best_per_run.csv, or, for a summary written before 04_summarize.py had
    them, the same rows rebuilt from all_k_per_run.csv.gz by the same rule
    (highest ``overall``, ties to the smaller K, as compare_ks.selection._best_k)."""
    rows = df[df["role"] == "annotation_best"]
    if len(rows):
        return rows
    logging.info("no annotation_best rows in best_per_run.csv - rebuilding them "
                 "from all_k_per_run.csv.gz")
    allk = pd.read_csv(paths.EXP1_SUMMARY / "all_k_per_run.csv.gz",
                       keep_default_na=False, na_values=[""])
    allk = allk[allk["ann_key"].notna()].sort_values("k")
    best = allk.loc[allk.groupby("run_id")["overall"].idxmax()]
    return (best.rename(columns={"sc": "reference", "overall": "score_overall"})
                .assign(role="annotation_best")[df.columns])


def best_start(df: pd.DataFrame) -> pd.DataFrame:
    """best_per_run.csv with each SpatialAIM row replaced by the better start,
    the SpatialAIM side of figures 2 and 4.

    Input
        df      best_per_run.csv: one row per run and role

    Output
        the baseline rows unchanged, plus one ``spatialaim`` row per (pair, tree,
        mapper): among the Leiden start and every annotation
        start, each at its own K = best.overall, the one with the highest
        ``score_overall``. Its ``ann_key`` says which start won (empty = Leiden).
    """
    join = ["pair_dir", "tree", "mapper"]
    candidates = pd.concat([df[df["role"] == "spatialaim"], annotation_best(df)],
                           ignore_index=True)
    best = candidates.loc[candidates.groupby(join)["score_overall"].idxmax()]
    won = best["ann_key"].isna()
    logging.info("best start: Leiden in %d of %d runs (%.0f%%)",
                 won.sum(), len(best), 100 * won.mean())
    best = best.assign(role="spatialaim")
    return pd.concat([df[df["role"] == "baseline"], best], ignore_index=True)


def paired(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    """Put each SpatialAIM run next to the annotation run it is compared with.

    Input
        df      best_start(best_per_run.csv): one row per run, role "spatialaim"
                (the better start, K = best.overall) or "baseline" (annotation
                start, K = T)
        metric  the metric to compare, e.g. "cossim_hard_raw_spot"

    Output
        one row per (pair, tree, mapper), i.e. per SpatialAIM run, with
          all columns of the SpatialAIM row (``metric``, ``k``, ...)
          <metric>_base, k_base, ann_key_base   the same for the baseline
          delta                                 metric - metric_base

    Both sides use the same pair, tree and mapper, so the difference is only the
    starting partition and K. Where a reference has two annotations, the baseline is
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
