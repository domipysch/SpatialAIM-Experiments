"""Put the per-K jobs of 05a_run_k.py together into the normal run folder, so that
03_status.py and 04_summarize.py see an ordinary finished run.

    python 05b_merge_k.py <tree>

Needs every K = 1..L done under $EXP1_ROOT/split_k/<run_id>/; otherwise it prints
the --array list of the missing ones and stops. It then

  * scores the K levels together with spatialaim's own cross-K comparison
    (k_comparison.csv, k_selection.json -- as a full sweep would have),
  * packs them into k_levels.zip, as run_one.py does,
  * writes resources.json for the whole run: wall and CPU time summed over the K
    jobs, peak memory the largest of them, and ``per_k`` with each job's own
    numbers (wall_s, peak_rss_tree_mb, seconds_map, seconds_analysis),
  * and publishes it as $EXP1_ROOT/<run_id>/ with DONE.

run.log only notes that the run was split; the per-K times are in resources.json's
``per_k``, and each K's own log stays in split_k/.
"""

from __future__ import annotations

import argparse
import datetime as dt
import importlib.util
import json
import os
import re
import shutil
import sys
import tempfile
from pathlib import Path

from spatialaim.analysis.compare_ks import compare_k_runs

import grid
import run_one

HERE = Path(__file__).resolve().parent
# 05a's name starts with a digit, so it is loaded by path
_spec = importlib.util.spec_from_file_location("run_k", HERE / "05a_run_k.py")
run_k = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(run_k)

LOG_LINE = re.compile(r"^(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3}) - \S+ - \w+ - (.*)$")
ALIGNER_CALL = re.compile(r"^ReferenceMapper\[\w+\] K=\d+$")
MAPPED = re.compile(r"^K=\s*\d+ mapped -> ")
SWEEP_END = re.compile(r"^Comparing K-runs")


def phase_times(log: Path) -> dict:
    """Mapping (aligner call -> mapped) and analysis (mapped -> comparison) of a
    one-K job, from its log's timestamps."""
    t = {}
    for line in log.read_text(encoding="utf-8", errors="replace").splitlines():
        m = LOG_LINE.match(line)
        if not m:
            continue
        stamp = dt.datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S,%f").timestamp()
        for name, pattern in (("call", ALIGNER_CALL), ("mapped", MAPPED), ("end", SWEEP_END)):
            if pattern.match(m.group(2)):
                t.setdefault(name, stamp)
    return {"seconds_map": round(t["mapped"] - t["call"], 3) if {"call", "mapped"} <= t.keys() else None,
            "seconds_analysis": round(t["end"] - t["mapped"], 3) if {"mapped", "end"} <= t.keys() else None}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tree", choices=grid.TREES)
    ap.add_argument("--unit", type=int, default=run_k.UNIT)
    ap.add_argument("--mapper", default=run_k.MAPPER, choices=grid.MAPPERS)
    args = ap.parse_args()

    unit = grid.units()[args.unit]
    row = grid.mapping_run(unit, args.tree, args.mapper)
    dest = run_one.ROOT / row["run_id"]
    if (dest / "DONE").exists():
        print(f"{row['run_id']} is already done")
        return 0
    split = run_k.split_dir(row)
    ks = list(range(1, run_k.n_levels(unit) + 1))
    missing = [k for k in ks if not (split / f"k_{k:03d}" / "DONE").exists()]
    if missing:
        print(f"{len(missing)} of {len(ks)} K not done yet, resubmit:\n"
              f"  sbatch.tinyfat --array={grid.array_spec(missing)} 05a_run_k.tinyfat.sbatch "
              f"{grid.root_arg(run_one.ROOT)} {args.tree}")
        return 1

    per_k = [json.loads((split / f"k_{k:03d}" / "resources.json").read_text(encoding="utf-8"))
             for k in ks]
    merged = Path(tempfile.mkdtemp(prefix="exp1merge_", dir=os.environ.get("TMPDIR")))
    staging = dest.parent / f".{dest.name}.tmp-{os.getpid()}"
    try:
        for k in ks:
            shutil.copytree(split / f"k_{k:03d}" / f"k_{k:03d}", merged / f"k_{k:03d}")
        compare_k_runs(merged, ks)

        if staging.exists():
            shutil.rmtree(staging)
        staging.mkdir(parents=True)
        shutil.copy2(split / f"k_{ks[0]:03d}" / "config.yaml", staging / "config.yaml")
        for name in ("k_comparison.csv", "k_selection.json"):
            shutil.copy2(merged / name, staging / name)
        zipped = run_one.zip_k_levels(merged, staging / run_one.K_ZIP)

        def total(key):
            values = [r.get(key) for r in per_k]
            return None if None in values else round(sum(values), 3)

        def peak(key):
            values = [r.get(key) for r in per_k if r.get(key) is not None]
            return max(values) if values else None

        res = {key: per_k[0].get(key) for key in (
            "run_id", "stage", "pair_dir", "start_tag", "ann_key", "tree", "mapper",
            "cluster", "partition", "cpu_model", "gpu_name", "threads", "mem_limit_mb",
            "spatialaim_version", "sample_interval_s")}
        res |= {
            "split_by_k": True,
            "started_at": min(r["started_at"] for r in per_k),
            "array_job_id": ";".join(sorted({str(r.get("array_job_id", "")) for r in per_k})),
            "host": ";".join(sorted({r.get("host", "") for r in per_k})),
            "exit_code": 0,
            "wall_s": total("wall_s"),
            "cpu_user_s": total("cpu_user_s"),
            "cpu_sys_s": total("cpu_sys_s"),
            "peak_rss_tree_mb": peak("peak_rss_tree_mb"),
            "peak_rss_proc_mb": peak("peak_rss_proc_mb"),
            "peak_rss_children_mb": peak("peak_rss_children_mb"),
            "peak_gpu_mem_mb": peak("peak_gpu_mem_mb"),
            "rss_samples": total("rss_samples"),
            "n_k": len(ks),
            "cache_hit": all(r.get("cache_hit") for r in per_k),
            "status": "ok",
            "per_k": [{"k": r["k"], "wall_s": r["wall_s"], "peak_rss_tree_mb": r["peak_rss_tree_mb"],
                       "host": r.get("host"), "job_id": r.get("job_id")}
                      | phase_times(split / f"k_{r['k']:03d}" / "run.log") for r in per_k],
        } | zipped
        (staging / "run.log").write_text(
            f"split into {len(ks)} one-K jobs (05a_run_k.py); their logs are in "
            f"{split.relative_to(run_one.ROOT).as_posix()}/k_<kkk>/run.log\n", encoding="utf-8")
        (staging / "resources.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
        (staging / "DONE").write_text(res["started_at"] + "\n")
        run_one._publish(staging, dest)
        run_one._ledger(res)
    finally:
        shutil.rmtree(merged, ignore_errors=True)
        shutil.rmtree(staging, ignore_errors=True)

    sel = json.loads((dest / "k_selection.json").read_text(encoding="utf-8"))
    print(f"merged {len(ks)} K -> {dest}\n  best overall K = {sel['best']['overall']}\n"
          f"  wall {res['wall_s'] / 3600:.1f} h summed over the K jobs")
    return 0


if __name__ == "__main__":
    sys.exit(main())
