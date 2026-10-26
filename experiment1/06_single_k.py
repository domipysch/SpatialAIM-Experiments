"""Runtime of mapping one annotation: per annotation unit and mapper, ``spatialaim run
--k_min T --k_max T`` on the unit's prebuilt scaffold, T = the number of annotation
cell types (the identity cut). Figure 2 sets it against a whole SpatialAIM run
(Leiden scaffold + every K).

Only the time is wanted; the results are those of the stage-B sweeps, which already
contain K = T. So a run keeps nothing but

    $EXP1_ROOT/single_k/<pair_dir>/<start_tag>/<mapper>/   run.log, resources.json, DONE | FAILED

and 04_summarize.py reads only its resources.json (-> summary/resources_single_k.csv).
It is measured like any run (run_one.measured, scaffold copied outside the timer),
so wall_s compares with stage B's; the K's own analysis runs, as it does for every
K of a sweep. One tree: at K = T the cut is every leaf, whatever the linkage.

    python 06_single_k.py <unit> <mapper>           # one array task (06_single_k.tinyfat.sbatch)
    python 06_single_k.py --missing                 # per mapper, the --array of the units without DONE
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

import grid
import run_one

HERE = Path(__file__).resolve().parent
TREE = "ward"


# 05a's name starts with a digit, so it is loaded by path
_spec = importlib.util.spec_from_file_location("run_k", HERE / "05a_run_k.py")
run_k = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(run_k)


def single_run(unit: dict, mapper: str) -> dict:
    return grid.mapping_run(unit, TREE, mapper) | {
        "run_id": f"single_k/{unit['pair_dir']}/{unit['start_tag']}/{mapper}"}


def run(row: dict) -> str:
    dest = run_one.ROOT / row["run_id"]
    if (dest / "DONE").exists():
        return "skip"
    scaffold = run_one.ROOT / grid.scaffold_run(row)["run_id"]
    if not (scaffold / "DONE").exists():
        return "no_scaffold"
    k = run_k.n_levels(row)

    scratch = Path(tempfile.mkdtemp(prefix="exp1t_", dir=os.environ.get("TMPDIR")))
    staging = dest.parent / f".{dest.name}.tmp-{os.getpid()}"
    try:
        out = scratch / "out"
        out.mkdir()
        for name in run_one.SCAFFOLD_FILES:  # outside the timer, as in run_one
            shutil.copy2(scaffold / name, out / name)
        cmd = run_one.command(row, out) + ["--k_min", str(k), "--k_max", str(k)]
        log_path = scratch / "run.log"
        res = run_one._context(row, cmd) | run_one.measured(cmd, log_path)
        log = log_path.read_text(encoding="utf-8", errors="replace")

        k_dir = out / row["mapper"] / f"k_{k:03d}"
        ok = all((k_dir / "analysis" / "data" / f).exists() for f in run_k.ANALYSIS_FILES)
        res |= {"k": k, "cache_hit": run_one.HIT in log, "status": "ok" if ok else "failed"}
        res |= run_one.phase_times(log_path)

        if staging.exists():
            shutil.rmtree(staging)
        staging.mkdir(parents=True)
        shutil.copy2(log_path, staging / "run.log")
        (staging / "resources.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
        (staging / ("DONE" if ok else "FAILED")).write_text(res["started_at"] + "\n")
        run_one._publish(staging, dest)
        run_one._ledger(res)
        return res["status"]
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
        shutil.rmtree(staging, ignore_errors=True)


def missing() -> None:
    annotation = [u for u in grid.units() if grid.kind(u) == "annotation"]
    for mapper in grid.MAPPERS:
        todo = [u["unit"] for u in annotation
                if not (run_one.ROOT / single_run(u, mapper)["run_id"] / "DONE").exists()]
        print(f"{mapper:<17} {len(annotation) - len(todo):>3} / {len(annotation)} done"
              + (f"   --array={grid.array_spec(todo)}%20" if todo else ""))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("unit", nargs="?", type=int, help="line of units.tsv (the array task id)")
    ap.add_argument("mapper", nargs="?", choices=grid.MAPPERS)
    ap.add_argument("--missing", action="store_true", help="what is left, per mapper")
    args = ap.parse_args()
    if args.missing:
        missing()
        return 0
    if args.unit is None or args.mapper is None:
        ap.error("give a unit and a mapper, or --missing")

    unit = grid.units()[args.unit]
    if grid.kind(unit) != "annotation":
        print(f"unit {args.unit} is a Leiden start; nothing to do")
        return 0
    row = single_run(unit, args.mapper)
    status = run(row)
    print(f"{status:<11} {row['run_id']}", flush=True)
    return 0 if status in ("ok", "skip") else 1


if __name__ == "__main__":
    sys.exit(main())
