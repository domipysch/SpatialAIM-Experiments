"""Run one demo-results job through the ``spatialaim`` CLI and publish its output
straight into the layout the online demo reads.

    scaffold <pair_id>      build the pair's reference scaffold -> both linkages' pair_dir
    sweep    <run>          full K sweep of one run (row of runs.tsv)
    k        <run> <K>      one K of a run split by K -> $DEMO_STATE/split_k/
    merge    [<run>]        put those K together into the run's normal folder (all split runs)
    levels   <run>          print the --array range 1-L of ``k`` for a run
    status                  done / failed / missing, with resubmit lines

Every job works in node-local scratch ($TMPDIR) with a copy of the pair's scaffold, so
concurrent jobs of one pair never write the same file; a finished run folder is moved
into $DEMO_OUT in one rename. Files shared by a pair's runs (``start_clustering.h5ad``,
``cell_shape_features.csv``) are written once, by the first run that has them. Markers
and logs go to $DEMO_STATE, so $DEMO_OUT holds results only: ``DONE`` = skip,
``FAILED`` or none = rerun on resubmit.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import time
import traceback
from pathlib import Path

import grid

OUT = Path(os.environ.get("DEMO_OUT", "results")).resolve()
STATE = Path(os.environ.get("DEMO_STATE", "state")).resolve()
DATA = Path(os.environ.get("DEMO_DATA", grid.data_dir())).resolve()
LEIDEN_RES = "3.0"

SCAFFOLD_FILES = ("reference_scaffold.h5ad", "reference_scaffold.meta.json")
SHAPE_FEATURES = "cell_shape_features.csv"  # spatialaim.analysis.shape_consistency.SHAPE_FEATURES_CSV
START_CLUSTERING = "start_clustering.h5ad"  # spatialaim.io.START_CLUSTERING_FILENAME
# the files the cross-K comparison reads; a K without them is not done
ANALYSIS_FILES = ("cossim_metrics.json", "topology_metrics.json", "modularity_metrics.json")


def pair_dir(row: dict, linkage: str | None = None) -> Path:
    return OUT / (linkage or row["linkage"]) / row["pair_dir"]


def run_dest(row: dict) -> Path:
    return pair_dir(row) / row["start_tag"] / row["mapper"]


def run_state(row: dict) -> Path:
    return STATE / "runs" / row["linkage"] / row["pair_dir"] / row["start_tag"] / row["mapper"]


def scaffold_state(row: dict) -> Path:
    return STATE / "scaffolds" / row["pair_dir"]


def split_dir(row: dict) -> Path:
    return STATE / "split_k" / row["linkage"] / row["pair_dir"] / row["start_tag"] / row["mapper"]


def _spatialaim() -> list[str]:
    exe = shutil.which("spatialaim")
    return [exe] if exe else [sys.executable, "-m", "spatialaim"]


def command(row: dict, out: Path, scaffold: bool = False, k: int | None = None) -> list[str]:
    cmd = _spatialaim() + [
        "run",
        "--scdata", str(DATA / "scRNA" / f"{row['sc']}.h5ad"),
        "--stdata", str(DATA / "ST" / f"{row['st']}.h5ad"),
        "--output_dir", str(out),
        "--mapping", "nearest_centroid" if scaffold else row["mapper"],
        "--linkage_method", "average" if scaffold else row["linkage"],
        "--leiden_resolution", LEIDEN_RES,
    ]
    if scaffold:
        return cmd + ["--k_min", "1", "--k_max", "1"]
    if row["ann_key"]:
        cmd += ["--annotation", row["ann_key"], "--no_leiden_start"]
    if row["segmentation"]:
        cmd += ["--segmentation", str(DATA / "ST" / "Segmentation" / row["segmentation"])]
    if k is not None:
        cmd += ["--k_min", str(k), "--k_max", str(k)]
    return cmd


def _copy_once(src: Path, dest: Path) -> None:
    """Copy ``src`` to ``dest`` unless it exists; the rename makes a concurrent
    reader see either nothing or the whole file."""
    if dest.exists() or not src.exists():
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(f".{dest.name}.tmp-{os.getpid()}")
    shutil.copy2(src, tmp)
    os.replace(tmp, dest)


def _publish_dir(src: Path, dest: Path) -> None:
    """Replace ``dest`` by a copy of ``src``, moved in with one rename."""
    staging = dest.with_name(f".{dest.name}.tmp-{os.getpid()}")
    shutil.rmtree(staging, ignore_errors=True)
    shutil.copytree(src, staging)
    if dest.exists():
        shutil.rmtree(dest)
    os.replace(staging, dest)


def _mark(state: Path, ok: bool, log: Path | None, note: str) -> None:
    state.mkdir(parents=True, exist_ok=True)
    for name in ("DONE", "FAILED"):
        (state / name).unlink(missing_ok=True)
    if log is not None and log.exists():
        shutil.copy2(log, state / "run.log")
    (state / ("DONE" if ok else "FAILED")).write_text(note + "\n", encoding="utf-8")


def _execute(cmd: list[str], scratch: Path) -> tuple[int, Path, float]:
    log = scratch / "run.log"
    t0 = time.time()
    with open(log, "wb") as fh:
        code = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT).returncode
    return code, log, time.time() - t0


def _scratch() -> Path:
    return Path(tempfile.mkdtemp(prefix="demo_", dir=os.environ.get("TMPDIR")))


def _with_scaffold(row: dict, out: Path) -> None:
    out.mkdir(parents=True)
    for name in SCAFFOLD_FILES:
        shutil.copy2(pair_dir(row, grid.LINKAGES[0]) / name, out / name)


def _publish_shared(row: dict, out: Path, linkages: list[str]) -> None:
    """The files a pair's runs share, into each linkage's pair_dir once."""
    for linkage in linkages:
        _copy_once(out / SHAPE_FEATURES, pair_dir(row, linkage) / SHAPE_FEATURES)
        if row.get("start_tag"):
            _copy_once(out / row["start_tag"] / START_CLUSTERING,
                       pair_dir(row, linkage) / row["start_tag"] / START_CLUSTERING)


def scaffold(row: dict) -> str:
    """Build the pair's scaffold (a K = 1 nearest_centroid run) and put it into every
    linkage's pair_dir."""
    state = scaffold_state(row)
    if (state / "DONE").exists():
        return "skip"
    scratch = _scratch()
    try:
        out = scratch / "out"
        code, log, secs = _execute(command(row, out, scaffold=True), scratch)
        ok = code == 0 and all((out / n).exists() for n in SCAFFOLD_FILES)
        if ok:
            for linkage in grid.LINKAGES:
                for name in SCAFFOLD_FILES:
                    _copy_once(out / name, pair_dir(row, linkage) / name)
        _mark(state, ok, log, f"exit {code}, {secs:.0f} s")
        return "ok" if ok else "failed"
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def sweep(row: dict) -> str:
    """The run's full K sweep into $DEMO_OUT."""
    state = run_state(row)
    if (state / "DONE").exists():
        return "skip"
    if not (scaffold_state(row) / "DONE").exists():
        return "no_scaffold"
    scratch = _scratch()
    try:
        out = scratch / "out"
        _with_scaffold(row, out)
        code, log, secs = _execute(command(row, out), scratch)
        root = out / row["start_tag"] / row["mapper"]
        ok = code == 0 and (root / "k_comparison.csv").exists()
        if ok:
            _publish_shared(row, out, [row["linkage"]])
            _publish_dir(root, run_dest(row))
        _mark(state, ok, log, f"exit {code}, {secs:.0f} s")
        return "ok" if ok else "failed"
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def n_levels(row: dict) -> int:
    """L, the largest K of a run's sweep: its number of start clusters."""
    import anndata as ad

    obs = ad.read_h5ad(pair_dir(row, grid.LINKAGES[0]) / SCAFFOLD_FILES[0], backed="r").obs
    return int(obs[row["ann_key"] or "leiden_all_genes"].nunique())


def one_k(row: dict, k: int) -> str:
    """One K of a run split by K, into split_k/<run>/k_<kkk>/."""
    dest = split_dir(row) / f"k_{k:03d}"
    marker = split_dir(row) / "done" / f"{k:03d}"
    if marker.exists() or (run_state(row) / "DONE").exists():
        return "skip"
    if not (scaffold_state(row) / "DONE").exists():
        return "no_scaffold"
    scratch = _scratch()
    try:
        out = scratch / "out"
        _with_scaffold(row, out)
        code, log, secs = _execute(command(row, out, k=k), scratch)
        root = out / row["start_tag"] / row["mapper"]
        k_dir = root / f"k_{k:03d}"
        ok = all((k_dir / "analysis" / f).exists() for f in ANALYSIS_FILES)
        log_dir = split_dir(row) / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(log, log_dir / f"k_{k:03d}.log")
        if ok:
            _publish_shared(row, out, [row["linkage"]])
            _copy_once(root / "config.yaml", split_dir(row) / "config.yaml")
            _publish_dir(k_dir, dest)
            marker.parent.mkdir(parents=True, exist_ok=True)
            marker.write_text(f"{secs:.0f} s\n", encoding="utf-8")
        return "ok" if ok else "failed"
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def merge(row: dict) -> int:
    """The K jobs of a split run as one ordinary run folder, compared across K by
    spatialaim's own K comparison."""
    from spatialaim.analysis.compare_ks import compare_k_runs

    if (run_state(row) / "DONE").exists():
        print(f"run {row['run']} is already done")
        return 0
    split = split_dir(row)
    ks = list(range(1, n_levels(row) + 1))
    missing = [k for k in ks if not (split / "done" / f"{k:03d}").exists()]
    if missing:
        print(f"{len(missing)} of {len(ks)} K not done yet, resubmit:\n"
              f"  sbatch.tinyfat --array={grid.array_spec(missing)} 02b_sweep_by_k.tinyfat.sbatch {row['run']}")
        return 1
    merged = Path(tempfile.mkdtemp(prefix="demomerge_", dir=os.environ.get("TMPDIR")))
    try:
        for k in ks:
            shutil.copytree(split / f"k_{k:03d}", merged / f"k_{k:03d}")
        # The config of a full sweep: no K range.
        config = (split / "config.yaml").read_text(encoding="utf-8").splitlines()
        config = [f"{line.split(':')[0]}: null" if line.startswith(("k_min:", "k_max:")) else line
                  for line in config]
        (merged / "config.yaml").write_text("\n".join(config) + "\n", encoding="utf-8")
        compare_k_runs(merged, ks)
        _publish_dir(merged, run_dest(row))
        _mark(run_state(row), True, None, f"merged from {len(ks)} one-K jobs")
    finally:
        shutil.rmtree(merged, ignore_errors=True)
    print(f"merged {len(ks)} K -> {run_dest(row)}")
    return 0


def status() -> int:
    """Done / failed / missing per job kind, and the lines to resubmit the rest."""
    rows = grid.runs()
    pairs = {r["pair_id"]: r for r in rows}
    todo_scaffold = [p for p, r in pairs.items() if not (scaffold_state(r) / "DONE").exists()]
    print(f"scaffolds: {len(pairs) - len(todo_scaffold)}/{len(pairs)} done")
    if todo_scaffold:
        print(f"  sbatch.tinyfat --array={grid.array_spec(todo_scaffold)} 01_scaffold.tinyfat.sbatch")
    open_runs = [r for r in rows if not (run_state(r) / "DONE").exists()]
    failed = [r for r in open_runs if (run_state(r) / "FAILED").exists()]
    print(f"runs: {len(rows) - len(open_runs)}/{len(rows)} done, {len(failed)} failed")
    for r in failed:
        print(f"  FAILED {r['run']:>4}  {run_state(r).relative_to(STATE)}  (log: run.log there)")
    for cls in grid.CLASSES:
        todo = [r["run"] for r in open_runs if grid.run_class(r) == cls]
        if todo:
            print(f"  sbatch.tinyfat -J demo_{cls} --time={grid.CLASSES[cls]} "
                  f"--array={grid.array_spec(todo)}%20 02_sweep.tinyfat.sbatch")
    for r in open_runs:
        if grid.split_by_k(r):
            print(f"  split by K, run {r['run']}: sbatch.tinyfat --array=$(python run_one.py levels {r['run']})%20 "
                  f"02b_sweep_by_k.tinyfat.sbatch {r['run']}")
    if any(grid.split_by_k(r) for r in open_runs):
        print("  once their K are done: python run_one.py merge")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("scaffold")
    p.add_argument("pair_id", type=int)
    p.add_argument("--dry-run", action="store_true", help="print the command only")
    p = sub.add_parser("sweep")
    p.add_argument("run", type=int, help="row of runs.tsv (the array task id)")
    p.add_argument("--dry-run", action="store_true", help="print the command only")
    p = sub.add_parser("k")
    p.add_argument("run", type=int)
    p.add_argument("k", type=int, help="the K (the array task id)")
    sub.add_parser("merge").add_argument("run", type=int, nargs="?",
                                         help="without: every run split by K")
    sub.add_parser("levels").add_argument("run", type=int)
    sub.add_parser("status")
    args = ap.parse_args()

    if args.cmd == "status":
        return status()
    if args.cmd == "merge" and args.run is None:
        return max((merge(r) for r in grid.runs() if grid.split_by_k(r)), default=0)
    if args.cmd == "scaffold":
        matches = [r for r in grid.runs() if int(r["pair_id"]) == args.pair_id]
        if not matches:
            ap.error(f"PairID {args.pair_id} is not in runs.tsv")
        row = matches[0]
        if args.dry_run:
            print(" ".join(command(row, Path("<scratch>/out"), scaffold=True)))
            return 0
        job = lambda: scaffold(row)  # noqa: E731
    else:
        row = grid.runs()[args.run]
        if args.cmd == "levels":
            print(f"1-{n_levels(row)}")
            return 0
        if args.cmd == "merge":
            return merge(row)
        if getattr(args, "dry_run", False):
            print(" ".join(command(row, Path("<scratch>/out"))))
            return 0
        job = (lambda: one_k(row, args.k)) if args.cmd == "k" else (lambda: sweep(row))

    t0 = time.time()
    try:
        result = job()
    except Exception:  # noqa: BLE001 - reported as a failed job
        traceback.print_exc()
        result = "error"
    print(f"{result:<11} {time.time() - t0:8.0f}s  {args.cmd} {vars(args)}", flush=True)
    return 0 if result in ("ok", "skip") else 1


if __name__ == "__main__":
    sys.exit(main())
