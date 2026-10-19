"""One K of one stage-B run, for the run too long for one 24 h job: DOT on unit 189
(pair 93, 09_human-breast-cancer, whole 168k-cell Xenium section, Leiden start).

The sweep is split into one Slurm array task per K (05a_run_k.tinyfat.sbatch):
each runs ``spatialaim run --k_min K --k_max K`` on the unit's prebuilt scaffold,
measured like any run (run_one.measured), and keeps that K's small files under

    $EXP1_ROOT/split_k/<run_id>/k_<kkk>/   k_<kkk>/..., config.yaml, run.log, resources.json, DONE

05b_merge_k.py then puts the K levels together into the normal run folder. A K
with DONE is skipped, so resubmitting the array only redoes what is missing.

    python 05a_run_k.py <tree> <K>
    python 05a_run_k.py <tree> --levels      # print the --array range 1-L

``--unit``/``--mapper`` exist for testing; the defaults are the case this is for.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

import anndata as ad

import grid
import run_one

UNIT, MAPPER = 189, "dot"
# the files the cross-K comparison reads; a K without them is not done
ANALYSIS_FILES = ("cossim_summary.csv", "topology_metrics.json", "modularity_metrics.json")


def split_dir(row: dict) -> Path:
    return run_one.ROOT / "split_k" / row["run_id"]


def n_levels(unit: dict) -> int:
    """L: the number of start clusters, i.e. the largest K."""
    scaffold = run_one.ROOT / grid.scaffold_run(unit)["run_id"] / run_one.SCAFFOLD_FILES[0]
    return int(ad.read_h5ad(scaffold, backed="r").uns["n_start_clusters"])


def keep_k(k_dir: Path, dest: Path) -> None:
    """Copy a K folder without the files that grow with the number of spots."""
    for path in sorted(k_dir.rglob("*")):
        if not path.is_file() or any(path.match(p) for p in run_one.PER_SPOT):
            continue
        if path.stat().st_size > run_one.ZIP_MAX_FILE_MB * 2**20:
            continue
        target = dest / path.relative_to(k_dir.parent)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)


def run_k(row: dict, k: int) -> str:
    dest = split_dir(row) / f"k_{k:03d}"
    if (dest / "DONE").exists():
        return "skip"
    scaffold = run_one.ROOT / grid.scaffold_run(row)["run_id"]
    if not (scaffold / "DONE").exists():
        return "no_scaffold"

    scratch = Path(tempfile.mkdtemp(prefix="exp1k_", dir=os.environ.get("TMPDIR")))
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

        root = out / row["mapper"]
        k_dir = root / f"k_{k:03d}"
        # A one-level sweep's own K selection is meaningless and may not even exit
        # cleanly; what counts is that the K was mapped and analysed.
        ok = all((k_dir / "analysis" / "data" / f).exists() for f in ANALYSIS_FILES)
        res |= {"k": k, "cache_hit": run_one.HIT in log, "status": "ok" if ok else "failed"}

        if staging.exists():
            shutil.rmtree(staging)
        staging.mkdir(parents=True)
        if ok:
            keep_k(k_dir, staging)
            shutil.copy2(root / "config.yaml", staging / "config.yaml")
        shutil.copy2(log_path, staging / "run.log")
        (staging / "resources.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
        (staging / ("DONE" if ok else "FAILED")).write_text(res["started_at"] + "\n")
        run_one._publish(staging, dest)
        return res["status"]
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
        shutil.rmtree(staging, ignore_errors=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tree", choices=grid.TREES)
    ap.add_argument("k", nargs="?", type=int, help="the K of this array task")
    ap.add_argument("--levels", action="store_true", help="print the --array range 1-L")
    ap.add_argument("--unit", type=int, default=UNIT)
    ap.add_argument("--mapper", default=MAPPER, choices=grid.MAPPERS)
    args = ap.parse_args()

    unit = grid.units()[args.unit]
    row = grid.mapping_run(unit, args.tree, args.mapper)
    if args.levels:
        print(f"1-{n_levels(unit)}")
        return 0
    if args.k is None:
        ap.error("give a K, or --levels")
    if (run_one.ROOT / row["run_id"] / "DONE").exists():
        print(f"{row['run_id']} is already done; nothing to split")
        return 0
    status = run_k(row, args.k)
    print(f"{status:<11} K={args.k}  {row['run_id']}", flush=True)
    return 0 if status in ("ok", "skip") else 1


if __name__ == "__main__":
    sys.exit(main())
