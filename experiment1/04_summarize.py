"""
Summarize all results from Experiment1 into a few csvs.
``summary/`` is what gets pulled.
This script only reads the run folders and edits/writes nothing outside ``summary/``.

Produces:
    summary/best_per_run.csv          per sweep the compared K: annotation start K = T and best overall, Leiden start best overall
    summary/all_k_per_run.csv.gz      every k_comparison.csv row, with the run's keys
    summary/resources_scaffold.csv    time + memory per scaffold build (one per pair), with its phases
    summary/resources_sweep.csv       time + memory per sweep, its phases, and the two derived costs
                                      full_pipeline_s / one_annotation_s (see sweep_costs)
    summary/resources_per_k.csv       per sweep and K: seconds and peak memory of the mapping and of the analysis
    summary/missing.csv               grid runs without DONE (should be empty)
    summary/meta.json                 counts, spatialaim version(s), when

"""

from __future__ import annotations
import datetime as dt
import json
import os
from pathlib import Path
import numpy as np
import pandas as pd
import grid

ROOT = Path(os.environ["EXP1_ROOT"])
OUT = ROOT / "summary"

# the keys every per-run table carries
KEYS = ["run_id", "pair_id", "pair_dir", "sc", "st", "start_tag", "ann_key", "tree", "mapper"]

# the metrics the figures plot: reconstruction (section 02), biology and topology (section 04),
# plus shape consistency (NaN for slices without a segmentation mask)
METRICS = ["cossim_hard_raw_spot", "cossim_hard_raw_gene",
           "spatial_assort_celltype", "spatial_assort_global",
           "modularity_st_expression", "modularity_shared", "shape"]

# the fields of resources.json that go into every resources table; the phase
# split (scaffold_s, setup_s, map_s, analysis_s, sweep_s, rest_s) is added as it comes
MEASURED = ["status", "wall_s", "cpu_user_s", "cpu_sys_s", "peak_rss_tree_mb",
            "peak_rss_proc_mb", "peak_rss_children_mb", "peak_gpu_mem_mb", "gpu_mem_base_mb",
            "rss_samples", "threads", "mem_limit_mb", "cluster", "partition", "host",
            "cpu_model", "gpu_name", "job_id", "spatialaim_version",
            "started_at", "cache_hit", "split_by_k", "n_k"]


def figure_rows(run: pd.Series) -> list[dict]:
    """The K of a run that the figures compare:

        baseline          annotation start, K = T (the identity cut)
        annotation_best   annotation start, K = best.overall of k_selection.json
        spatialaim        Leiden start, K = best.overall of k_selection.json
    """
    table = pd.read_csv(ROOT / run["run_id"] / "k_comparison.csv")
    table = table.set_index(table["k"].astype(int))
    selection = json.loads((ROOT / run["run_id"] / "k_selection.json").read_text(encoding="utf-8"))
    best = (selection.get("best") or {}).get("overall")
    best_k = int(best["k"]) if best is not None and int(best["k"]) in table.index else None
    if run["ann_key"]:
        roles = [("baseline", int(table.index.max())), ("annotation_best", best_k)]
    else:
        roles = [("spatialaim", best_k)]
    return [{"pair_id": run["pair_id"], "reference": run["sc"], "pair_dir": run["pair_dir"],
             "role": role, "ann_key": run["ann_key"], "tree": run["tree"], "mapper": run["mapper"],
             "k": k, "score_overall": float(table.loc[k, "overall"]),
             **{m: float(table.loc[k].get(m, float("nan"))) for m in METRICS}}
            for role, k in roles if k is not None]


def measured(run: pd.Series) -> dict:
    """``resources.json`` of a run: MEASURED plus every phase total (``*_s``)."""
    try:
        res = json.loads((ROOT / run["run_id"] / "resources.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        res = {}
    phases = {f: v for f, v in res.items() if f.endswith("_s") and f not in MEASURED
              and f != "sample_interval_s"}
    return {f: res.get(f) for f in MEASURED} | phases | (
        {"scaffold_mb": res.get("scaffold_mb")} if run["stage"] == "scaffold" else {})


def per_k(run: pd.Series) -> pd.DataFrame:
    """One row per K of a sweep, from its ``phases.csv``: seconds and peak memory
    of ``map`` (the mapper call) and ``analysis`` (writing + analysing that K)."""
    phases = pd.read_csv(ROOT / run["run_id"] / "phases.csv", dtype={"k": "Int64"})
    phases = phases[phases["k"].notna() & phases["phase"].isin(["map", "analysis"])]
    wide = phases.pivot_table(index="k", columns="phase",
                              values=["seconds", "peak_rss_tree_mb", "peak_rss_children_mb"],
                              aggfunc="max")
    wide.columns = [f"{value}_{phase}" for value, phase in wide.columns]
    return wide.reset_index().assign(**run[KEYS])


def sweep_costs(sweeps: pd.DataFrame, scaffolds: pd.DataFrame, k_rows: pd.DataFrame) -> pd.DataFrame:
    """Two derived costs per sweep, from measured times only.

    scaffold_build_s    the pair's scaffold build (``scaffold_s`` of its scaffold job)
    full_pipeline_s     the sweep as if it had built its own scaffold: its wall time
                        with its scaffold phase (the cache read) replaced by the build
    one_annotation_s    annotation starts only: mapping that one annotation at K = T
                        and analysing it -- the sweep's wall time without the K below
                        T (it sweeps from T down, so T is its first K)
    """
    build = scaffolds.set_index("pair_dir")["scaffold_s"].reindex(sweeps["pair_dir"]).to_numpy()
    out = sweeps.assign(scaffold_build_s=build)
    out["full_pipeline_s"] = out["wall_s"] - out["scaffold_s"] + build
    at_t = (k_rows.sort_values("k").groupby("run_id").tail(1).set_index("run_id")
            [["seconds_map", "seconds_analysis"]].sum(axis=1).reindex(out["run_id"]).to_numpy())
    one = out["wall_s"] - out["sweep_s"] + at_t
    out["one_annotation_s"] = np.where(out["ann_key"].astype(bool), one, np.nan)
    return out


def main() -> None:

    # Load
    OUT.mkdir(parents=True, exist_ok=True)
    runs = pd.DataFrame(grid.runs(ROOT))
    runs["pair_id"] = runs["pair_id"].astype(int)
    runs["done"] = [(ROOT / r / "DONE").exists() for r in runs["run_id"]]
    sweeping = runs[(runs["stage"] == "sweep") & runs["done"]]
    building = runs[(runs["stage"] == "scaffold") & runs["done"]]

    # missing.csv
    runs.loc[~runs["done"], KEYS + ["stage"]].to_csv(OUT / "missing.csv", index=False)

    # best_per_run.csv
    rows = [r for _, run in sweeping.iterrows() for r in figure_rows(run)]
    figure = (pd.DataFrame(rows)
              .sort_values(["pair_id", "role", "ann_key", "tree", "mapper"], na_position="first")
              .reset_index(drop=True))
    figure.to_csv(OUT / "best_per_run.csv", index=False)

    # all_k_per_run.csv
    curves = [pd.read_csv(ROOT / run["run_id"] / "k_comparison.csv").assign(**run[KEYS])
              for _, run in sweeping.iterrows()]
    pd.concat(curves, ignore_index=True).to_csv(OUT / "all_k_per_run.csv.gz", index=False)

    # resources_scaffold.csv
    pair_keys = ["run_id", "pair_id", "pair_dir", "sc", "st"]
    scaffolds = pd.DataFrame([dict(run[pair_keys]) | measured(run) for _, run in building.iterrows()])
    scaffolds.to_csv(OUT / "resources_scaffold.csv", index=False)

    # resources_per_k.csv
    k_rows = pd.concat([per_k(run) for _, run in sweeping.iterrows()], ignore_index=True)
    k_rows = k_rows[KEYS + [c for c in k_rows.columns if c not in KEYS]]
    k_rows.to_csv(OUT / "resources_per_k.csv", index=False)

    # resources_sweep.csv
    sweeps = pd.DataFrame([dict(run[KEYS]) | measured(run) for _, run in sweeping.iterrows()])
    sweeps = sweep_costs(sweeps, scaffolds, k_rows)
    sweeps.to_csv(OUT / "resources_sweep.csv", index=False)

    # meta.json
    both = pd.concat([scaffolds, sweeps])
    info = {
        "created": dt.datetime.now().isoformat(timespec="seconds"),
        "root": str(ROOT),
        "grid_runs": int(len(runs)),
        "done": int(runs["done"].sum()),
        "missing": int((~runs["done"]).sum()),
        "cache_misses": int((sweeps["cache_hit"].astype(str).str.lower() != "true").sum()),
        "spatialaim_versions": sorted(both["spatialaim_version"].dropna().unique().tolist()),
    }
    (OUT / "meta.json").write_text(json.dumps(info, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
