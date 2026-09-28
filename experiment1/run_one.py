"""Run one Experiment 1 job through the ``spatialaim`` CLI, measuring each run:
wall time, CPU time, peak memory of the whole process tree (the ``conda run``
aligner children included), and peak GPU memory where a GPU is visible (none
on TinyFat; kept for a run on a GPU node).

A job is one unit (line of ``units.tsv``, see ``grid.py``) in one stage:

    A <unit>            build the reference scaffold  -> scaffolds/<pair_dir>/<start_tag>
    B <unit> <mapper>   full K sweep, ward then average
                                          -> results/<pair_dir>/<start_tag>/<tree>/<mapper>

Stage A runs ``nearest_centroid`` at K = 1 only and keeps just
``reference_scaffold.{h5ad,meta.json}``: the scaffold is mapper- and
linkage-independent, so all ten stage-B runs of that (pair, start) reuse it. Each
stage-B run gets its own fresh output dir with a copy of the scaffold, so its
measurement is the sweep without the scaffold build; the log is checked for
the cache hit.

Every run works in node-local scratch (``$TMPDIR``) and copies back only the
small files: stage A the scaffold; stage B the run-root tables and one
``k_levels.zip`` with every per-K file except those that grow with the number
of spots (the mappings and the per-spot tables, see ``PER_SPOT``). A finished run has ``DONE`` and is skipped; a failed one has
``FAILED`` and is retried on the next call.

    python run_one.py A 0
    python run_one.py B 0 wann [--tree ward] [--dry-run]
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import platform
import re
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
SAMPLE_S = 0.5
GPU_SAMPLE_S = 1.0

HIT = "Reusing cached reference scaffold"
SCAFFOLD_WRITTEN = "Wrote reference scaffold cache"
SCAFFOLD_FILES = ("reference_scaffold.h5ad", "reference_scaffold.meta.json")
RUN_ROOT_FILES = ("config.yaml", "k_comparison.csv", "k_selection.json")
K_ZIP = "k_levels.zip"
# Per-K files left out of the zip: they scale with the number of spots.
PER_SPOT = ("*.h5ad", "spot_to_state_mapping*.csv", "onehot_per_row_mapping.csv",
            "cossim-per-spot-*", "confidence_per_spot.csv")
ZIP_MAX_FILE_MB = 2  # any other per-K file above this is left out too (and counted)
LOG_TIME = re.compile(r"^(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}),(\d{3}) - ")


# ---------------------------------------------------------------------------
# measurement
# ---------------------------------------------------------------------------

class Monitor:
    """Samples the RSS of a process and all its descendants until it exits."""

    def __init__(self, pid: int):
        self.root = psutil.Process(pid)
        self.peak_tree = 0
        self.peak_proc = 0
        self.peak_children = 0
        self.samples = 0
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
                self.peak_tree = max(self.peak_tree, sum(rss))
                self.peak_proc = max(self.peak_proc, max(rss))
                self.peak_children = max(self.peak_children, sum(rss[1:]))
                self.samples += 1
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


def measured(cmd: list[str], log_path: Path) -> dict:
    """Run ``cmd`` with stdout+stderr to ``log_path``; return the measurements."""
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
    mb = 1024 * 1024
    return {
        "exit_code": code,
        "started_at": dt.datetime.fromtimestamp(started).isoformat(timespec="seconds"),
        "started_epoch": started,
        "wall_s": round(wall, 3),
        "cpu_user_s": None if cpu0 is None else round(cpu1[0] - cpu0[0], 3),
        "cpu_sys_s": None if cpu0 is None else round(cpu1[1] - cpu0[1], 3),
        "peak_rss_tree_mb": round(mon.peak_tree / mb, 1),
        "peak_rss_proc_mb": round(mon.peak_proc / mb, 1),
        "peak_rss_children_mb": round(mon.peak_children / mb, 1),
        "peak_gpu_mem_mb": mon.peak_gpu,
        "gpu_mem_base_mb": mon.gpu_base,
        "rss_samples": mon.samples,
        "sample_interval_s": SAMPLE_S,
    }


def _log_epoch(line: str) -> float | None:
    m = LOG_TIME.match(line)
    if not m:
        return None
    stamp = dt.datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S")
    return stamp.timestamp() + int(m.group(2)) / 1000


# ---------------------------------------------------------------------------
# one run
# ---------------------------------------------------------------------------

def _spatialaim() -> list[str]:
    exe = shutil.which("spatialaim")
    return [exe] if exe else [sys.executable, "-m", "spatialaim"]


def command(row: dict, out: Path) -> list[str]:
    stage_a = row["stage"] == "A"
    cmd = _spatialaim() + [
        "run",
        "--scdata", str(DATA / "scRNA" / f"{row['sc']}.h5ad"),
        "--stdata", str(DATA / "ST" / f"{row['st']}.h5ad"),
        "--output_dir", str(out),
        "--mapping", "nearest_centroid" if stage_a else row["mapper"],
        "--linkage_method", "average" if stage_a else row["tree"],
        "--leiden_resolution", LEIDEN_RES,
    ]
    if row["ann_key"]:
        cmd += ["--start_from_annotation", row["ann_key"]]
    if stage_a:
        cmd += ["--k_min", "1", "--k_max", "1"]
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
    """Pack the small files of every ``k_*/`` under ``root`` into ``dest``
    (paths kept relative to ``root``, e.g. ``k_012/analysis/data/biology_metrics.json``)."""
    n_files = n_skipped_large = 0
    with zipfile.ZipFile(dest, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(root.glob("k_*/**/*")):
            if not path.is_file() or any(path.match(p) for p in PER_SPOT):
                continue
            if path.stat().st_size > ZIP_MAX_FILE_MB * 2**20:
                n_skipped_large += 1
                continue
            zf.write(path, path.relative_to(root).as_posix())
            n_files += 1
    return {"zip_files": n_files, "zip_skipped_large": n_skipped_large,
            "zip_mb": round(dest.stat().st_size / 2**20, 2)}


def _publish(staging: Path, dest: Path) -> None:
    """Replace ``dest`` by ``staging`` (same filesystem, so near-atomic)."""
    if dest.exists():
        shutil.rmtree(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    os.replace(staging, dest)


def run(row: dict) -> str:
    dest = ROOT / row["run_id"]
    if (dest / "DONE").exists():
        return "skip"
    stage_a = row["stage"] == "A"
    scaffold = ROOT / "scaffolds" / row["pair_dir"] / row["start_tag"]
    if not stage_a and not (scaffold / "DONE").exists():
        return "no_scaffold"

    scratch = Path(tempfile.mkdtemp(prefix="exp1_", dir=os.environ.get("TMPDIR")))
    staging = dest.parent / f".{dest.name}.tmp-{os.getpid()}"
    try:
        out = scratch / "out"
        out.mkdir()
        if not stage_a:  # outside the timer: the run starts with its scaffold in place
            for name in SCAFFOLD_FILES:
                shutil.copy2(scaffold / name, out / name)
        cmd = command(row, out)
        log_path = scratch / "run.log"
        res = _context(row, cmd) | measured(cmd, log_path)
        log = log_path.read_text(encoding="utf-8", errors="replace")

        keep = []
        if stage_a:
            ok = all((out / n).exists() for n in SCAFFOLD_FILES)
            keep = [(out / n, n) for n in SCAFFOLD_FILES]
            written = next((_log_epoch(l) for l in log.splitlines() if SCAFFOLD_WRITTEN in l), None)
            res["scaffold_s"] = None if written is None else round(written - res["started_epoch"], 3)
            res["scaffold_mb"] = round((out / SCAFFOLD_FILES[0]).stat().st_size / 2**20, 1) if ok else None
        else:
            root = out / row["mapper"]
            ok = res["exit_code"] == 0 and (root / "k_comparison.csv").exists()
            keep = [(root / n, n) for n in RUN_ROOT_FILES if (root / n).exists()]
            if ok:
                zipped = zip_k_levels(root, scratch / K_ZIP)
                keep.append((scratch / K_ZIP, K_ZIP))
                res |= zipped
            res["cache_hit"] = HIT in log
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
        (staging / "resources.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
        (staging / ("DONE" if ok else "FAILED")).write_text(res["started_at"] + "\n")
        _publish(staging, dest)
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

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=["A", "B"])
    ap.add_argument("unit", type=int, help="line of units.tsv (the array task id)")
    ap.add_argument("mapper", nargs="?", choices=grid.MAPPERS, help="stage B only")
    ap.add_argument("--tree", choices=grid.TREES, help="stage B: only this linkage")
    ap.add_argument("--dry-run", action="store_true", help="print the commands only")
    args = ap.parse_args()
    if (args.stage == "B") != (args.mapper is not None):
        ap.error("a mapper is required for stage B and not allowed for stage A")

    unit = grid.units()[args.unit]
    if args.stage == "A":
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
