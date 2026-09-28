"""
Summarize all results from Experiment1 into a few csvs.
``summary/`` is what gets pulled.
This script only reads the run folders and edits/writes nothing outside ``summary/``.

Produces:
    summary/best_per_run.csv          one row per stage-B run, (start-from-annotation: K=T, overclustering: best overall score)
    summary/all_k_per_run.csv.gz      every k_comparison.csv row, with the run's keys
    summary/resources_stage_b.csv     time + memory per stage-B run, with its scaffold's cost
    summary/resources_stage_a.csv     time + memory per stage-A scaffold build
    summary/missing.csv               grid runs without DONE (should be empty)
    summary/meta.json                 counts, spatialaim version(s), when

"""

from __future__ import annotations
import datetime as dt
import json
import os
from pathlib import Path
import pandas as pd
import grid

ROOT = Path(os.environ["EXP1_ROOT"])
OUT = ROOT / "summary"

# the keys every per-run table carries
KEYS = ["run_id", "pair_id", "pair_dir", "sc", "st", "start_tag", "ann_key", "tree", "mapper"]

# the metrics the figures plot: reconstruction (section 03), biology and topology (section 04)
METRICS = ["cossim_hard_raw_spot", "cossim_hard_raw_gene",
           "nhood_mean_self_zscore", "local_purity_zscore",
           "modularity_st_expression", "modularity_shared"]

# the fields of resources.json that go into resources_stage_b.csv / resources_stage_a.csv
MEASURED = ["status", "wall_s", "cpu_user_s", "cpu_sys_s", "peak_rss_tree_mb",
            "peak_rss_proc_mb", "peak_rss_children_mb", "peak_gpu_mem_mb", "gpu_mem_base_mb",
            "rss_samples", "threads", "mem_limit_mb", "cluster", "partition", "host",
            "cpu_model", "gpu_name", "job_id", "spatialaim_version", "started_at"]
SCAFFOLD_ONLY = ["scaffold_s", "scaffold_mb"]


def figure_row(run: pd.Series) -> dict | None:
    """The one K of a run that the figures compare:

        baseline     annotation start, K = T (the identity cut)
        spatialaim   Leiden start, K = best.overall of k_selection.json
    """
    table = pd.read_csv(ROOT / run["run_id"] / "k_comparison.csv")
    table = table.set_index(table["k"].astype(int))
    if run["ann_key"]:
        role, k = "baseline", int(table.index.max())
    else:
        selection = json.loads((ROOT / run["run_id"] / "k_selection.json").read_text(encoding="utf-8"))
        best = (selection.get("best") or {}).get("overall")
        if best is None or int(best["k"]) not in table.index:
            return None
        role, k = "spatialaim", int(best["k"])
    row = table.loc[k]
    return {"pair_id": run["pair_id"], "reference": run["sc"], "pair_dir": run["pair_dir"],
            "role": role, "ann_key": run["ann_key"], "tree": run["tree"], "mapper": run["mapper"],
            "k": k, "score_overall": float(row["overall"]),
            **{m: float(row.get(m, float("nan"))) for m in METRICS}}


def measured(run: pd.Series, fields: list[str]) -> dict:
    try:
        res = json.loads((ROOT / run["run_id"] / "resources.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        res = {}
    return {f: res.get(f) for f in fields}


def main() -> None:

    # Load
    OUT.mkdir(parents=True, exist_ok=True)
    runs = pd.DataFrame(grid.runs(ROOT))
    runs["pair_id"] = runs["pair_id"].astype(int)
    runs["done"] = [(ROOT / r / "DONE").exists() for r in runs["run_id"]]
    mapping = runs[(runs["stage"] == "B") & runs["done"]]
    building = runs[(runs["stage"] == "A") & runs["done"]]

    # missing.csv
    runs.loc[~runs["done"], KEYS + ["stage"]].to_csv(OUT / "missing.csv", index=False)

    # best_per_run.csv
    figure_rows = [r for _, run in mapping.iterrows() if (r := figure_row(run)) is not None]
    figure = (pd.DataFrame(figure_rows)
              .sort_values(["pair_id", "role", "ann_key", "tree", "mapper"], na_position="first")
              .reset_index(drop=True))
    figure.to_csv(OUT / "best_per_run.csv", index=False)

    # all_k_per_run.csv
    curves = [pd.read_csv(ROOT / run["run_id"] / "k_comparison.csv").assign(**run[KEYS])
              for _, run in mapping.iterrows()]
    pd.concat(curves, ignore_index=True).to_csv(OUT / "all_k_per_run.csv.gz", index=False)

    # resources_stage_a.csv
    scaffold_keys = [k for k in KEYS if k not in ("tree", "mapper")]
    scaffolds = pd.DataFrame([dict(run[scaffold_keys]) | measured(run, MEASURED + SCAFFOLD_ONLY)
                              for _, run in building.iterrows()])
    scaffolds.to_csv(OUT / "resources_stage_a.csv", index=False)

    # resources_stage_b.csv
    # each mapping run with the cost of the scaffold it reused: the build share of
    # stage A (log-timed), its whole wall time where that is missing
    cost = scaffolds.set_index(["pair_dir", "start_tag"])
    resources = pd.DataFrame([dict(run[KEYS]) | measured(run, MEASURED)
                              for _, run in mapping.iterrows()])
    unit = pd.MultiIndex.from_frame(resources[["pair_dir", "start_tag"]])
    resources["scaffold_build_s"] = cost["scaffold_s"].fillna(cost["wall_s"]).reindex(unit).to_numpy()
    resources["scaffold_peak_rss_tree_mb"] = cost["peak_rss_tree_mb"].reindex(unit).to_numpy()
    resources["wall_with_scaffold_s"] = resources["wall_s"] + resources["scaffold_build_s"]
    resources.to_csv(OUT / "resources_stage_b.csv", index=False)

    # meta.json
    versions = pd.concat([scaffolds["spatialaim_version"], resources["spatialaim_version"]])
    info = {
        "created": dt.datetime.now().isoformat(timespec="seconds"),
        "root": str(ROOT),
        "grid_runs": int(len(runs)),
        "done": int(runs["done"].sum()),
        "missing": int((~runs["done"]).sum()),
        "spatialaim_versions": sorted(versions.dropna().unique().tolist()),
    }
    (OUT / "meta.json").write_text(json.dumps(info, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
