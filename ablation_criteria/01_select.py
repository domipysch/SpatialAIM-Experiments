"""
The K each criterion set selects on every Leiden-start sweep of Experiment 1.

The sets are those of ``k_selection.json``: each selection criterion alone, every
combination of 2 to n - 1 of them, and all of them. The selection is SpatialAIM's
own (``compare_ks.best_ks`` / ``best_combination_ks``) on the sweep's scores in
``all_k_per_run.csv.gz``, so "all" is the K of ``best_per_run.csv``.

Writes ``<ABLATION_CRITERIA>/selected_k.csv``: one row per sweep and set with the
sweep's keys, ``criteria``, ``k``, ``score`` and ``k_max`` (the sweep's largest K).
"""

from __future__ import annotations
import logging
import pandas as pd
from spatialaim.analysis.compare_ks import best_combination_ks, best_ks
from spatialaim.analysis.registry import SELECTION_CRITERIA
from spatialaim_experiments import paths

KEYS = ["run_id", "pair_id", "pair_dir", "sc", "st", "start_tag", "tree", "mapper"]


def selected(table: pd.DataFrame) -> dict[str, tuple[int, float] | None]:
    """The best ``(k, score)`` of one sweep per criterion set."""
    table = table.sort_values("k", ignore_index=True)
    return best_ks(table, SELECTION_CRITERIA) | best_combination_ks(table, SELECTION_CRITERIA)


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    allk = pd.read_csv(paths.EXP1_SUMMARY / "all_k_per_run.csv.gz",
                       keep_default_na=False, na_values=[""])
    allk = allk[allk["start_tag"] == "leiden"]
    rows = []
    for _, run in allk.groupby("run_id", sort=False):
        keys = run.iloc[0][KEYS].to_dict() | {"k_max": int(run["k"].max())}
        for criteria, hit in selected(run).items():
            rows.append(keys | {"criteria": criteria,
                                "k": hit[0] if hit else None,
                                "score": hit[1] if hit else None})
    out = pd.DataFrame(rows).astype({"k": "Int64"})

    # "all" must be the K Experiment 1 reports
    best = pd.read_csv(paths.EXP1_SUMMARY / "best_per_run.csv")
    best = best[best["role"] == "spatialaim"].set_index(["pair_dir", "tree", "mapper"])["k"]
    mine = out[out["criteria"] == "all"].set_index(["pair_dir", "tree", "mapper"])["k"]
    mismatch = (mine.reindex(best.index) != best).sum()
    assert mismatch == 0, f"{mismatch} sweeps select another K than best_per_run.csv"

    paths.ABLATION_CRITERIA.mkdir(parents=True, exist_ok=True)
    out.to_csv(paths.ABLATION_CRITERIA / "selected_k.csv", index=False)
    logging.info("%d sweeps x %d sets -> %s", out["run_id"].nunique(),
                 out["criteria"].nunique(), paths.ABLATION_CRITERIA / "selected_k.csv")


if __name__ == "__main__":
    main()
