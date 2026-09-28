"""The Experiment 1 grid, fixed up front in ``units.tsv`` (checked in).

A *unit* is one (pair, start clustering); line i of ``units.tsv`` is unit i and
also the Slurm array index of every job that works on it:

    stage A   job (unit)           builds the reference scaffold       01_scaffolds.tinyfat.sbatch
    stage B   job (unit, mapper)   runs ward, then average             02_map.tinyfat.sbatch <mapper>

Start clusterings per pair: ``leiden_res3.0`` plus ``ann_<key>`` for every
``CellTypeKey*`` of the reference -> 252 units for the 88 pairs. Units are
ordered by PairID, so units 0-2 are pair 1 (01_mouse-ssp); the test runs units 0
(Leiden start) and 1 (ann_cellType), see TEST_UNITS.

    python grid.py array leiden|annotation|all [--test]   # the --array list of one start kind
    python grid.py units                                  # rewrite units.tsv from $EXP1_DATA (only if the datasets change)

Run folders, relative to the run root:

    scaffolds/<pair_dir>/<start_tag>
    results/<pair_dir>/<start_tag>/<tree>/<mapper>
"""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
UNITS = HERE / "units.tsv"
LEIDEN_TAG = "leiden_res3.0"
TREES = ["ward", "average"]
MAPPERS = ["nearest_centroid", "wann", "tacco", "dot", "tangram"]
# The one-pair test: pair 1 (01_mouse-ssp), one Leiden and one annotation start.
TEST_UNITS = ("0", "1")

# Walltime per stage-B job (ward + average), from the first round's slowest jobs:
# nearest_centroid, wann, tangram (CPU) < 1 h, tacco 1.3 h -> 2 h, which also
# qualifies for TinyFat's short-job nodes. DOT depends on the start kind
# (annotation 2.6 h p99 / 5.1 h max, leiden 6.6 / 8.2 h); its limits are realistic
# rather than worst-case, and the rare job that runs over is resubmitted with
# --time=24:00:00 (its finished tree has DONE and is skipped).
JOB_TIME = "02:00:00"
DOT_TIME = {"annotation": "04:00:00", "leiden": "08:00:00"}


def units() -> list[dict]:
    with open(UNITS, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def kind(unit: dict) -> str:
    return "leiden" if unit["start_tag"] == LEIDEN_TAG else "annotation"


def job_time(mapper: str, unit: dict) -> str:
    return DOT_TIME[kind(unit)] if mapper == "dot" else JOB_TIME


def scaffold_run(unit: dict) -> dict:
    """The stage-A run of a unit."""
    return unit | {"run_id": f"scaffolds/{unit['pair_dir']}/{unit['start_tag']}",
                   "stage": "A", "tree": "", "mapper": ""}


def mapping_run(unit: dict, tree: str, mapper: str) -> dict:
    """One stage-B run of a unit."""
    return unit | {"run_id": f"results/{unit['pair_dir']}/{unit['start_tag']}/{tree}/{mapper}",
                   "stage": "B", "tree": tree, "mapper": mapper}


def runs(root: Path) -> list[dict]:
    """Every run a run root is meant to hold: per unit one scaffold and ten mapping
    runs; the test root (``.../test``) only TEST_UNITS."""
    chosen = [u for u in units() if root.name != "test" or u["unit"] in TEST_UNITS]
    return [run for u in chosen
            for run in [scaffold_run(u)] + [mapping_run(u, t, m) for m in MAPPERS for t in TREES]]


def array_spec(indices) -> str:
    """[0, 1, 2, 5, 7, 8] -> '0-2,5,7-8'"""
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
    """How the sbatch scripts name ``root``: run, test, or its absolute path."""
    return root.name if root.name in ("run", "test") and root.parent.name == "spatialaim_exp1" else str(root)


def write_units(data: Path) -> None:
    """units.tsv from pairs.csv and the references' CellTypeKey* columns."""
    with open(data / "scRNA" / "index.csv", newline="", encoding="utf-8") as fh:
        keys = {row["Name"]: [k for k in ((row.get(f"CellTypeKey{i}") or "").strip() for i in range(3))
                              if k and k.lower() != "nan"]
                for row in csv.DictReader(fh)}
    with open(data / "pairs.csv", newline="", encoding="utf-8") as fh:
        pairs = sorted(csv.DictReader(fh), key=lambda p: int(p["PairID"]))
    rows = []
    for p in pairs:
        pair_dir = f"{int(p['PairID']):03d}_{p['scName']}__{p['stName']}"
        for start_tag, ann_key in [(LEIDEN_TAG, "")] + [(f"ann_{k}", k) for k in keys[p["scName"]]]:
            rows.append({"unit": len(rows), "pair_id": int(p["PairID"]), "pair_dir": pair_dir,
                         "sc": p["scName"], "st": p["stName"],
                         "start_tag": start_tag, "ann_key": ann_key})
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
