"""Run one Experiment 1 job through the ``spatialaim`` CLI, measuring it: wall
time, CPU time, and the memory of the whole process tree (the ``conda run``
aligner children included), sampled every SAMPLE_S seconds and kept as a time
series. spatialaim's own ``timing.csv`` (the scaffold, and per K the mapping
and the analysis) is joined with that series, so each of them gets its seconds
and its peak memory (``phases.csv``); the rest of the run's time follows from
the process start and end (``phase_totals``).

    scaffold <unit>                     build the pair's reference scaffold -> scaffolds/<pair_dir>
    sweep    <unit> <mapper> [--tree]   full K sweep, ward then average
                                        -> results/<pair_dir>/<start_tag>/<tree>/<mapper>
    k        <unit> <mapper> <tree> <K> one K of a sweep too long for one job -> split_k/<run_id>/k_<kkk>
    merge    <unit> <mapper> <tree>     put those K together into the normal run folder
    levels   <unit>                     print the --array range 1-L of ``k`` for a unit

The scaffold job runs the Leiden start with ``nearest_centroid`` at K = 1 and
keeps only ``reference_scaffold.{h5ad,meta.json}``: the scaffold is start-,
mapper- and linkage-independent, so every sweep of the pair reuses it. Each
sweep gets its own fresh output dir with a copy of the scaffold (outside the
timer), so its measurement holds no scaffold build; ``cache_hit`` checks that
(spatialaim's log line on reusing the cache).

Every run works in node-local scratch (``$TMPDIR``) and copies back only the
small files: the scaffold job the scaffold; a sweep the run-root tables and one
``k_levels.zip`` with every per-K file except the mappings (``*.h5ad``, which
grow with the number of spots). Each run also keeps ``run.log``,
``resources.json``, ``phases.csv`` and ``memory.tsv.gz``. A finished run has
``DONE`` and is skipped; a failed one has ``FAILED`` and is retried.

    python run_one.py scaffold 0
    python run_one.py sweep 0 wann [--tree ward] [--dry-run]
"""

from __future__ import annotations

import argparse
import bisect
import csv
import datetime as dt
import gzip
import json
import os
import platform
import shutil
import socket
import subprocess
import sys
import tempfile
import threading
import time
import traceback
import zipfile
from pathlib import Path

import psutil

import grid

try:
    import resource  # POSIX only: CPU time of all waited-for descendants
except ImportError:  # pragma: no cover - Windows smoke tests
    resource = None

ROOT = Path(os.environ.get("EXP1_ROOT", ".")).resolve()
DATA = Path(os.environ.get("EXP1_DATA", ".")).resolve()
LEIDEN_RES = "3.0"
SAMPLE_S = 0.25
GPU_SAMPLE_S = 1.0
MB = 2**20

SCAFFOLD_FILES = ("reference_scaffold.h5ad", "reference_scaffold.meta.json")
RUN_ROOT_FILES = ("config.yaml", "k_comparison.csv", "k_selection.json")
TIMING = "timing.csv"  # spatialaim.timing.TIMING_CSV
K_ZIP = "k_levels.zip"
# the files the cross-K comparison reads; a K without them is not done
ANALYSIS_FILES = ("cossim_metrics.json", "topology_metrics.json", "modularity_metrics.json")
HIT = "Reusing cached reference scaffold"  # spatialaim.io.read_reference_scaffold


# ---------------------------------------------------------------------------
# measurement
# ---------------------------------------------------------------------------

class Monitor:
    """Samples the RSS of a process and all its descendants until it exits:
    the peaks, and the series ``(epoch, tree MB, children MB)``."""

    def __init__(self, pid: int):
        self.root = psutil.Process(pid)
        self.peak_tree = 0
        self.peak_proc = 0
        self.peak_children = 0
        self.series: list[tuple[float, float, float]] = []
        self.peak_gpu = None
        self.gpu_base = None
        self._stop = threading.Event()
        self._threads = [threading.Thread(target=self._rss, daemon=True)]
        if _gpu_memory() is not None:
            self.gpu_base = max(_gpu_memory())
            self.peak_gpu = self.gpu_base
            self._threads.append(threading.Thread(target=self._gpu, daemon=True))
        for t in self._threads:
            t.start()

    def _rss(self) -> None:
        while not self._stop.is_set():
            try:
                procs = [self.root] + self.root.children(recursive=True)
            except psutil.NoSuchProcess:
                break
            rss = []
            for p in procs:
                try:
                    rss.append(p.memory_info().rss)
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
            if rss:
                tree, children = sum(rss), sum(rss[1:])
                self.peak_tree = max(self.peak_tree, tree)
                self.peak_proc = max(self.peak_proc, max(rss))
                self.peak_children = max(self.peak_children, children)
                self.series.append((time.time(), tree / MB, children / MB))
            self._stop.wait(SAMPLE_S)

    def _gpu(self) -> None:
        while not self._stop.is_set():
            used = _gpu_memory()
            if used:
                self.peak_gpu = max(self.peak_gpu, max(used))
            self._stop.wait(GPU_SAMPLE_S)

    def stop(self) -> None:
        self._stop.set()
        for t in self._threads:
            t.join(timeout=5)


def _gpu_memory() -> list[int] | None:
    """MiB in use per visible GPU, or None without a GPU."""
    if not os.environ.get("CUDA_VISIBLE_DEVICES") or not shutil.which("nvidia-smi"):
        return None
    try:
        out = subprocess.run(
            ["nvidia-smi", "--query-gpu=memory.used", "--format=csv,noheader,nounits"],
            capture_output=True, text=True, timeout=10, check=True,
        ).stdout
        return [int(v) for v in out.split()]
    except Exception:  # noqa: BLE001 - a failed sample is not a failed run
        return None


def _gpu_name() -> str:
    if _gpu_memory() is None:
        return ""
    out = subprocess.run(["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
                         capture_output=True, text=True).stdout
    return ";".join(line.strip() for line in out.splitlines() if line.strip())


def _cpu_model() -> str:
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor()


def _mem_limit_mb() -> int | None:
    if os.environ.get("SLURM_MEM_PER_NODE"):
        return int(os.environ["SLURM_MEM_PER_NODE"])
    if os.environ.get("SLURM_MEM_PER_CPU") and os.environ.get("SLURM_CPUS_PER_TASK"):
        return int(os.environ["SLURM_MEM_PER_CPU"]) * int(os.environ["SLURM_CPUS_PER_TASK"])
    return None


def _children_cpu() -> tuple[float, float] | None:
    if resource is None:
        return None
    ru = resource.getrusage(resource.RUSAGE_CHILDREN)
    return ru.ru_utime, ru.ru_stime


def measured(cmd: list[str], log_path: Path) -> tuple[dict, list]:
    """Run ``cmd`` with stdout+stderr to ``log_path``; return the measurements and
    the memory series."""
    cpu0 = _children_cpu()
    started = time.time()
    t0 = time.perf_counter()
    with open(log_path, "wb") as log:
        proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT)
        mon = Monitor(proc.pid)
        code = proc.wait()
    wall = time.perf_counter() - t0
    mon.stop()
    cpu1 = _children_cpu()
    return {
        "exit_code": code,
        "started_at": dt.datetime.fromtimestamp(started).isoformat(timespec="seconds"),
        "started_epoch": started,
        "wall_s": round(wall, 3),
        "cpu_user_s": None if cpu0 is None else round(cpu1[0] - cpu0[0], 3),
        "cpu_sys_s": None if cpu0 is None else round(cpu1[1] - cpu0[1], 3),
        "peak_rss_tree_mb": round(mon.peak_tree / MB, 1),
        "peak_rss_proc_mb": round(mon.peak_proc / MB, 1),
        "peak_rss_children_mb": round(mon.peak_children / MB, 1),
        "peak_gpu_mem_mb": mon.peak_gpu,
        "gpu_mem_base_mb": mon.gpu_base,
        "rss_samples": len(mon.series),
        "sample_interval_s": SAMPLE_S,
    }, mon.series


def write_series(series: list, path: Path) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
        fh.write("epoch\trss_tree_mb\trss_children_mb\n")
        for t, tree, children in series:
            fh.write(f"{t:.2f}\t{tree:.1f}\t{children:.1f}\n")


def phases(timing_files: list[tuple[str, Path]], series: list) -> list[dict]:
    """spatialaim's ``timing.csv`` rows (``scope`` = which file), each with the
    peak memory of the samples inside its window. A window shorter than the
    sampling interval may hold no sample; it gets the first sample after its
    start (``peak_sampled`` False)."""
    epochs = [s[0] for s in series]
    rows = []
    for scope, path in timing_files:
        if not path.exists():
            continue
        with open(path, newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                start, seconds = float(row["start_epoch"]), float(row["seconds"])
                lo = bisect.bisect_left(epochs, start)
                hi = bisect.bisect_right(epochs, start + seconds)
                sampled = hi > lo
                inside = series[lo:hi] if sampled else series[lo:lo + 1]
                rows.append({
                    "scope": scope, "phase": row["phase"], "k": row["k"],
                    "start_epoch": round(start, 3), "seconds": round(seconds, 3),
                    "peak_rss_tree_mb": max((s[1] for s in inside), default=None),
                    "peak_rss_children_mb": max((s[2] for s in inside), default=None),
                    "peak_sampled": sampled,
                })
    return rows


def phase_totals(rows: list[dict], res: dict) -> dict:
    """The run's time split, ``setup_s + sweep_s + rest_s = wall_s``:

        setup_s    process start to the first K's mapping (imports, loading, shape
                   features, the scaffold -- also on its own as ``scaffold_s`` --,
                   start clusters, tree, mapper preparation)
        sweep_s    map + analysis over all K (also as ``map_s`` and ``analysis_s``)
        rest_s     after the last K: closing the mapper, the K comparison, exit
    """
    def total(phase):
        return sum(r["seconds"] for r in rows if r["phase"] == phase)

    maps = [r for r in rows if r["phase"] == "map"]
    first = min((r["start_epoch"] for r in maps), default=None)
    setup = None if first is None else first - res["started_epoch"]
    sweep = total("map") + total("analysis")
    return {
        "scaffold_s": round(total("scaffold"), 3),
        "setup_s": None if setup is None else round(setup, 3),
        "map_s": round(total("map"), 3),
        "analysis_s": round(total("analysis"), 3),
        "sweep_s": round(sweep, 3),
        "rest_s": None if setup is None else round(res["wall_s"] - setup - sweep, 3),
    }


def write_rows(rows: list[dict], path: Path) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, list(rows[0]) if rows else ["scope"])
        w.writeheader()
        w.writerows(rows)


# ---------------------------------------------------------------------------
# one run
# ---------------------------------------------------------------------------

def _spatialaim() -> list[str]:
    exe = shutil.which("spatialaim")
    return [exe] if exe else [sys.executable, "-m", "spatialaim"]


def command(row: dict, out: Path, k: int | None = None) -> list[str]:
    scaffold = row["stage"] == "scaffold"
    cmd = _spatialaim() + [
        "run",
        "--scdata", str(DATA / "scRNA" / f"{row['sc']}.h5ad"),
        "--stdata", str(DATA / "ST" / f"{row['st']}.h5ad"),
        "--output_dir", str(out),
        "--mapping", "nearest_centroid" if scaffold else row["mapper"],
        "--linkage_method", "average" if scaffold else row["tree"],
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


def _context(row: dict, cmd: list[str]) -> dict:
    try:
        from importlib.metadata import version
        spatialaim_version = version("spatial-aim")
    except Exception:  # noqa: BLE001
        spatialaim_version = ""
    return {
        "run_id": row["run_id"],
        "stage": row["stage"],
        "pair_dir": row["pair_dir"],
        "start_tag": row["start_tag"],
        "ann_key": row["ann_key"],
        "tree": row["tree"],
        "mapper": row["mapper"],
        "cluster": os.environ.get("SLURM_CLUSTER_NAME", ""),
        "partition": os.environ.get("SLURM_JOB_PARTITION", ""),
        "job_id": os.environ.get("SLURM_JOB_ID", ""),
        "array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID", ""),
        "array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID", ""),
        "host": socket.gethostname(),
        "cpu_model": _cpu_model(),
        "gpu_name": _gpu_name(),
        "threads": int(os.environ.get("OMP_NUM_THREADS", "0") or 0),
        "mem_limit_mb": _mem_limit_mb(),
        "spatialaim_version": spatialaim_version,
        "command": " ".join(cmd),
    }


def zip_k_levels(root: Path, dest: Path) -> dict:
    """Pack every ``k_*/`` file under ``root`` except the mappings (``*.h5ad``)
    into ``dest``, paths relative to ``root`` (``k_012/analysis/topology_metrics.json``)."""
    n_files = 0
    with zipfile.ZipFile(dest, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(root.glob("k_*/**/*")):
            if path.is_file() and path.suffix != ".h5ad":
                zf.write(path, path.relative_to(root).as_posix())
                n_files += 1
    return {"zip_files": n_files, "zip_mb": round(dest.stat().st_size / MB, 2)}


def _publish(staging: Path, dest: Path) -> None:
    """Replace ``dest`` by ``staging`` (same filesystem, so near-atomic)."""
    if dest.exists():
        shutil.rmtree(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    os.replace(staging, dest)


def run(row: dict, k: int | None = None) -> str:
    """One run (``row``, see grid.py), or with ``k`` only that K of a sweep,
    into ``split_k/<run_id>/k_<kkk>``."""
    dest = ROOT / row["run_id"] if k is None else split_dir(row) / f"k_{k:03d}"
    if (dest / "DONE").exists():
        return "skip"
    is_scaffold = row["stage"] == "scaffold"
    scaffold = ROOT / grid.scaffold_run(row)["run_id"]
    if not is_scaffold and not (scaffold / "DONE").exists():
        return "no_scaffold"

    scratch = Path(tempfile.mkdtemp(prefix="exp1_", dir=os.environ.get("TMPDIR")))
    staging = dest.parent / f".{dest.name}.tmp-{os.getpid()}"
    try:
        out = scratch / "out"
        out.mkdir()
        if not is_scaffold:  # outside the timer: the run starts with its scaffold in place
            for name in SCAFFOLD_FILES:
                shutil.copy2(scaffold / name, out / name)
        cmd = command(row, out, k)
        log_path = scratch / "run.log"
        res, series = measured(cmd, log_path)
        res = _context(row, cmd) | res
        res["cache_hit"] = HIT in log_path.read_text(encoding="utf-8", errors="replace")
        start_dir = out / (grid.LEIDEN_TAG if is_scaffold else row["start_tag"])
        root = start_dir / ("nearest_centroid" if is_scaffold else row["mapper"])
        rows = phases([("pair", out / TIMING), ("run", root / TIMING)], series)
        res |= phase_totals(rows, res)

        if is_scaffold:
            ok = res["exit_code"] == 0 and all((out / n).exists() for n in SCAFFOLD_FILES)
            keep = [(out / n, n) for n in SCAFFOLD_FILES]
            if ok:
                res["scaffold_mb"] = round((out / SCAFFOLD_FILES[0]).stat().st_size / MB, 1)
        elif k is not None:
            # A one-level sweep's own K selection is meaningless; what counts is
            # that the K was mapped and analysed.
            k_dir = root / f"k_{k:03d}"
            ok = all((k_dir / "analysis" / f).exists() for f in ANALYSIS_FILES)
            keep = [(p, p.relative_to(root).as_posix()) for p in sorted(k_dir.rglob("*"))
                    if p.is_file() and p.suffix != ".h5ad"]
            keep.append((root / "config.yaml", "config.yaml"))
            res["k"] = k
        else:
            ok = res["exit_code"] == 0 and (root / "k_comparison.csv").exists()
            keep = [(root / n, n) for n in RUN_ROOT_FILES if (root / n).exists()]
            if ok:
                res |= zip_k_levels(root, scratch / K_ZIP)
                keep.append((scratch / K_ZIP, K_ZIP))
            res["n_k"] = len(list(root.glob("k_*/")))
        res["status"] = "ok" if ok else "failed"

        if staging.exists():
            shutil.rmtree(staging)
        staging.mkdir(parents=True)
        if ok:
            for src, rel in keep:
                (staging / rel).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, staging / rel)
        shutil.copy2(log_path, staging / "run.log")
        write_rows(rows, staging / "phases.csv")
        write_series(series, staging / "memory.tsv.gz")
        (staging / "resources.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
        (staging / ("DONE" if ok else "FAILED")).write_text(res["started_at"] + "\n")
        _publish(staging, dest)
        if k is None:
            _ledger(res)
        return res["status"]
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
        shutil.rmtree(staging, ignore_errors=True)  # only left behind by an error


def _ledger(res: dict) -> None:
    fields = ["started_at", "status", "wall_s", "peak_rss_tree_mb", "peak_gpu_mem_mb",
              "exit_code", "host", "job_id", "run_id"]
    path = ROOT / "logs" / "ledger.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    line = "\t".join(str(res.get(f, "")) for f in fields) + "\n"
    with open(path, "a", encoding="utf-8") as fh:  # one short O_APPEND write per run
        fh.write(line)


# ---------------------------------------------------------------------------
# a sweep split by K (one job per K), for a run that does not fit one job
# ---------------------------------------------------------------------------

def split_dir(row: dict) -> Path:
    return ROOT / "split_k" / row["run_id"]


def n_levels(unit: dict) -> int:
    """L, the largest K of a unit's sweep: its number of start clusters."""
    import anndata as ad
    obs = ad.read_h5ad(ROOT / grid.scaffold_run(unit)["run_id"] / SCAFFOLD_FILES[0], backed="r").obs
    return int(obs[unit["ann_key"] or "leiden_all_genes"].nunique())


def merge(row: dict, unit: dict) -> int:
    """The K jobs of ``row`` as one ordinary run folder: the K levels compared
    together by spatialaim's own cross-K comparison, the times summed over the
    jobs, the peak memory the largest of them, and every job's phases in
    ``phases.csv`` (``job`` = its K)."""
    from spatialaim.analysis.compare_ks import compare_k_runs

    dest = ROOT / row["run_id"]
    if (dest / "DONE").exists():
        print(f"{row['run_id']} is already done")
        return 0
    split = split_dir(row)
    ks = list(range(1, n_levels(unit) + 1))
    missing = [k for k in ks if not (split / f"k_{k:03d}" / "DONE").exists()]
    if missing:
        print(f"{len(missing)} of {len(ks)} K not done yet, resubmit:\n"
              f"  sbatch.tinyfat --array={grid.array_spec(missing)} 02b_sweep_by_k.tinyfat.sbatch "
              f"{grid.root_arg(ROOT)} {unit['unit']} {row['mapper']} {row['tree']}")
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
        zipped = zip_k_levels(merged, staging / K_ZIP)
        rows = []
        for k in ks:
            with open(split / f"k_{k:03d}" / "phases.csv", newline="", encoding="utf-8") as fh:
                rows += [{"job": k} | r for r in csv.DictReader(fh)]
        write_rows(rows, staging / "phases.csv")

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
        summed = sorted({key for r in per_k for key in r
                         if key.endswith("_s") and key != "sample_interval_s"})
        res |= {
            "split_by_k": True,
            "started_at": min(r["started_at"] for r in per_k),
            "array_job_id": ";".join(sorted({str(r.get("array_job_id", "")) for r in per_k})),
            "host": ";".join(sorted({r.get("host", "") for r in per_k})),
            "exit_code": 0,
            **{key: total(key) for key in summed},
            "peak_rss_tree_mb": peak("peak_rss_tree_mb"),
            "peak_rss_proc_mb": peak("peak_rss_proc_mb"),
            "peak_rss_children_mb": peak("peak_rss_children_mb"),
            "peak_gpu_mem_mb": peak("peak_gpu_mem_mb"),
            "rss_samples": total("rss_samples"),
            "n_k": len(ks),
            "cache_hit": all(r.get("cache_hit") for r in per_k),
            "status": "ok",
        } | zipped
        (staging / "run.log").write_text(
            f"split into {len(ks)} one-K jobs (run_one.py k); their logs and memory series are in "
            f"{split.relative_to(ROOT).as_posix()}/k_<kkk>/\n", encoding="utf-8")
        (staging / "resources.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
        (staging / "DONE").write_text(res["started_at"] + "\n")
        _publish(staging, dest)
        _ledger(res)
    finally:
        shutil.rmtree(merged, ignore_errors=True)
        shutil.rmtree(staging, ignore_errors=True)

    sel = json.loads((dest / "k_selection.json").read_text(encoding="utf-8"))
    print(f"merged {len(ks)} K -> {dest}\n  best overall K = {sel['best']['overall']}\n"
          f"  wall {res['wall_s'] / 3600:.1f} h summed over the K jobs")
    return 0


# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("scaffold", help="build the reference scaffold of a unit's pair")
    p.add_argument("unit", type=int, help="line of units.tsv (the array task id)")
    p.add_argument("--dry-run", action="store_true", help="print the command only")
    p = sub.add_parser("sweep", help="full K sweep of a unit, ward then average")
    p.add_argument("unit", type=int)
    p.add_argument("mapper", choices=grid.MAPPERS)
    p.add_argument("--tree", choices=grid.TREES, help="only this linkage")
    p.add_argument("--dry-run", action="store_true", help="print the commands only")
    p = sub.add_parser("k", help="one K of a sweep (a run split by K)")
    p.add_argument("unit", type=int)
    p.add_argument("mapper", choices=grid.MAPPERS)
    p.add_argument("tree", choices=grid.TREES)
    p.add_argument("k", type=int, help="the K (the array task id)")
    p = sub.add_parser("merge", help="merge a run split by K")
    p.add_argument("unit", type=int)
    p.add_argument("mapper", choices=grid.MAPPERS)
    p.add_argument("tree", choices=grid.TREES)
    p = sub.add_parser("levels", help="print the --array range 1-L of `k` for a unit")
    p.add_argument("unit", type=int)
    args = ap.parse_args()

    unit = grid.units()[args.unit]
    if args.cmd == "levels":
        print(f"1-{n_levels(unit)}")
        return 0
    if args.cmd == "merge":
        return merge(grid.mapping_run(unit, args.tree, args.mapper), unit)
    if args.cmd == "k":
        row = grid.mapping_run(unit, args.tree, args.mapper)
        if (ROOT / row["run_id"] / "DONE").exists():
            print(f"{row['run_id']} is already done; nothing to split")
            return 0
        status = run(row, args.k)
        print(f"{status:<11} K={args.k}  {row['run_id']}", flush=True)
        return 0 if status in ("ok", "skip") else 1

    if args.cmd == "scaffold":
        if grid.kind(unit) != "leiden":
            ap.error("the scaffold job runs on a pair's Leiden unit (python grid.py array leiden)")
        rows = [grid.scaffold_run(unit)]
    else:
        rows = [grid.mapping_run(unit, tree, args.mapper)
                for tree in ([args.tree] if args.tree else grid.TREES)]

    failed = 0
    for n, row in enumerate(rows, 1):
        if args.dry_run:
            print(" ".join(command(row, Path("<scratch>/out"))))
            continue
        t0 = time.time()
        try:
            status = run(row)
        except Exception:  # noqa: BLE001 - one broken run must not end the job
            traceback.print_exc()
            status = "error"
        failed += status not in ("ok", "skip")
        print(f"[{n}/{len(rows)}] {status:<11} {time.time() - t0:8.0f}s  {row['run_id']}", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
