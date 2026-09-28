"""Progress of one Experiment 1 run root: done / failed / missing per stage and
mapper, cache misses, and for everything unfinished the ``--array=`` list to
resubmit (a resubmitted job skips its runs that already have DONE).

    python 03_status.py              # table + resubmit lists
    python 03_status.py --failed     # + one line per failed run (exit code, peak memory, log)

Reads ``units.tsv`` (see ``grid.py``) and the DONE/FAILED markers under
``$EXP1_ROOT``. Pure stdlib + pandas.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import pandas as pd

import grid


def state(root: Path, run_id: str) -> str:
    d = root / run_id
    if (d / "DONE").exists():
        return "done"
    if (d / "FAILED").exists():
        return "failed"
    return "missing"


def resources(root: Path, run_id: str) -> dict:
    try:
        return json.loads((root / run_id / "resources.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def cause(res: dict) -> str:
    limit, peak = res.get("mem_limit_mb"), res.get("peak_rss_tree_mb") or 0
    if res.get("exit_code") in (-9, 137) or (limit and peak > 0.9 * limit):
        return "OOM?"
    return f"exit {res.get('exit_code')}"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--failed", action="store_true", help="list every failed run")
    args = ap.parse_args()
    root = Path(os.environ["EXP1_ROOT"])

    df = pd.DataFrame(grid.runs(root))
    df["unit"] = df["unit"].astype(int)
    df["state"] = [state(root, r) for r in df["run_id"]]
    df["mapper"] = df["mapper"].where(df["stage"] == "B", "(scaffold)")
    table = (df.groupby(["stage", "mapper", "state"]).size()
             .unstack(fill_value=0).reindex(columns=["done", "failed", "missing"], fill_value=0))
    table["total"] = table.sum(axis=1)
    print(f"root {root}\n")
    print(table.to_string())
    print(f"\nall: {int((df['state'] == 'done').sum())} / {len(df)} done")

    done_b = df[(df["stage"] == "B") & (df["state"] == "done")]
    misses = [r for r in done_b["run_id"] if resources(root, r).get("cache_hit") is False]
    print(f"cache misses among done stage-B runs: {len(misses)}")
    for r in misses[:20]:
        print(f"  {r}")

    todo = df[df["state"] != "done"]
    if len(todo):
        name = grid.root_arg(root)
        print("\nunfinished units, to resubmit (from the experiment1 dir):")
        todo = todo.assign(time=[
            "" if r["stage"] == "A" else grid.job_time(r["mapper"], r)
            for r in todo.to_dict("records")])
        for (stage, mapper, time), part in todo.groupby(["stage", "mapper", "time"]):
            spec = grid.array_spec(part["unit"].unique())
            if stage == "A":
                print(f"  sbatch.tinyfat --array={spec} 01_scaffolds.tinyfat.sbatch {name}")
            else:
                print(f"  sbatch.tinyfat -J exp1_B_{mapper} --time={time} --array={spec} "
                      f"02_map.tinyfat.sbatch {name} {mapper}")
        print("  (a job that ran out of time: resubmit it with a longer --time before the script\n"
              "   name, e.g. --time=02:00:00 for a scaffold, --time=24:00:00 for dot)")

    failed = df[df["state"] == "failed"]
    if args.failed and len(failed):
        print("\nfailed:")
        for r in failed["run_id"]:
            res = resources(root, r)
            print(f"  {cause(res):<8} {res.get('peak_rss_tree_mb', '?'):>9} MB  "
                  f"{res.get('wall_s', '?'):>9} s  {root / r / 'run.log'}")
    elif len(failed):
        print(f"\n{len(failed)} failed - rerun with --failed to list them")


if __name__ == "__main__":
    main()
