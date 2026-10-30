"""The Experiment 1 grid: one unit per (pair, start), listed in ``units.tsv``.

Row i of ``units.tsv`` is unit i and the Slurm array index of every job on it.
Each pair gets one scaffold job (on its Leiden unit) and per unit ten sweeps
(2 trees x 5 mappers).

    python grid.py array leiden|annotation|all [--test]   # print the --array list of one start kind
    python grid.py units                                  # rewrite units.tsv from $EXP1_DATA

Run folders, relative to the run root:

    scaffolds/<pair_dir>
    results/<pair_dir>/<start_tag>/<tree>/<mapper>
"""

from __future__ import annotations

import argparse
import csv
import os
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
UNITS = HERE / "units.tsv"
LEIDEN_TAG = "leiden"
TREES = ["ward", "average"]
MAPPERS = ["nearest_centroid", "wann", "tacco", "dot", "tangram"]
# The one-pair test: pair 1, its Leiden start and its cellType start.
TEST_UNITS = ("0", "1")

# Walltime per sweep job, from the first round's slowest jobs; a job that runs over is resubmitted.
JOB_TIME = "02:00:00"
DOT_TIME = {"annotation": "04:00:00", "leiden": "08:00:00"}


def start_tag(ann_key: str) -> str:
    """The start folder name, as spatialaim.io.start_folder_name builds it."""
    return LEIDEN_TAG if not ann_key else "annotation_" + re.sub(r"[^A-Za-z0-9._-]+", "_", ann_key)


def units() -> list[dict]:
    """All rows of units.tsv."""
    with open(UNITS, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def kind(unit: dict) -> str:
    """The unit's start kind: "leiden" or "annotation"."""
    return "leiden" if unit["start_tag"] == LEIDEN_TAG else "annotation"


def job_time(mapper: str, unit: dict) -> str:
    """The walltime of a sweep job."""
    return DOT_TIME[kind(unit)] if mapper == "dot" else JOB_TIME


def scaffold_run(unit: dict) -> dict:
    """The scaffold build of the unit's pair."""
    return unit | {"run_id": f"scaffolds/{unit['pair_dir']}",
                   "stage": "scaffold", "tree": "", "mapper": ""}


def mapping_run(unit: dict, tree: str, mapper: str) -> dict:
    """One sweep of the unit."""
    return unit | {"run_id": f"results/{unit['pair_dir']}/{unit['start_tag']}/{tree}/{mapper}",
                   "stage": "sweep", "tree": tree, "mapper": mapper}


def runs(root: Path) -> list[dict]:
    """Every run a run root should hold (the test root: TEST_UNITS only)."""
    chosen = [u for u in units() if root.name != "test" or u["unit"] in TEST_UNITS]
    return ([scaffold_run(u) for u in chosen if kind(u) == "leiden"]
            + [mapping_run(u, t, m) for u in chosen for m in MAPPERS for t in TREES])


def array_spec(indices) -> str:
    """Indices as a Slurm --array list: [0, 1, 2, 5] -> '0-2,5'."""
    indices = sorted(int(i) for i in indices)
    parts, start = [], None
    for i, u in enumerate(indices):
        if start is None:
            start = u
        if i + 1 == len(indices) or indices[i + 1] != u + 1:
            parts.append(str(start) if start == u else f"{start}-{u}")
            start = None
    return ",".join(parts)


def root_arg(root: Path) -> str:
    """The run root as the sbatch scripts take it: run, test, or its absolute path."""
    return root.name if root.name in ("run", "test") and root.parent.name == "spatialaim_exp1" else str(root)


def write_units(data: Path) -> None:
    """Rewrite units.tsv from pairs.csv and the two index.csv files."""
    # each reference's annotations (CellTypeKey0-2)
    with open(data / "scRNA" / "index.csv", newline="", encoding="utf-8") as fh:
        keys = {row["Name"]: [k for k in ((row.get(f"CellTypeKey{i}") or "").strip() for i in range(3))
                              if k and k.lower() != "nan"]
                for row in csv.DictReader(fh)}
    # each slice's segmentation mask, if any
    with open(data / "ST" / "index.csv", newline="", encoding="utf-8") as fh:
        masks = {row["Name"]: (row.get("Segmentation") or "").strip() for row in csv.DictReader(fh)}
    with open(data / "pairs.csv", newline="", encoding="utf-8") as fh:
        pairs = sorted(csv.DictReader(fh), key=lambda p: int(p["PairID"]))
    # per pair: the Leiden unit, then one unit per annotation
    rows = []
    for p in pairs:
        pair_dir = f"{int(p['PairID']):03d}_{p['scName']}__{p['stName']}"
        for ann_key in [""] + keys[p["scName"]]:
            rows.append({"unit": len(rows), "pair_id": int(p["PairID"]), "pair_dir": pair_dir,
                         "sc": p["scName"], "st": p["stName"],
                         "start_tag": start_tag(ann_key), "ann_key": ann_key,
                         "segmentation": masks.get(p["stName"], "")})
    with open(UNITS, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} units from {len(pairs)} pairs -> {UNITS}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("units", help="rewrite units.tsv from $EXP1_DATA")
    arr = sub.add_parser("array", help="print the --array index list of one start kind")
    arr.add_argument("kind", choices=["leiden", "annotation", "all"])
    arr.add_argument("--test", action="store_true", help="TEST_UNITS only (the test)")
    args = ap.parse_args()
    if args.cmd == "units":
        write_units(Path(os.environ.get("EXP1_DATA", HERE.parents[1] / "Data" / "01_Datasets_Paper")))
        return
    print(array_spec(u["unit"] for u in units()
                     if args.kind in ("all", kind(u)) and (not args.test or u["unit"] in TEST_UNITS)))


if __name__ == "__main__":
    main()
