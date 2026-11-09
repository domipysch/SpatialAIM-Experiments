"""The demo-results grid: one run per (pair, start, linkage, mapper), listed in ``runs.tsv``.

Row i of ``runs.tsv`` is run i and the Slurm array index of its sweep job; a scaffold
job's array index is the PairID. Starts are Leiden plus every ``CellTypeKey*`` of the
reference, as in Experiment 1.

    python grid.py runs                 # rewrite runs.tsv from $DEMO_DATA
    python grid.py pairs                # --array list of the scaffold jobs (the PairIDs)
    python grid.py array <class>        # --array list of one walltime class
    python grid.py time <class>         # its walltime

Output, relative to $DEMO_OUT:  <linkage>/<pair_dir>/<start_tag>/<mapper>/
"""

from __future__ import annotations

import argparse
import csv
import os
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs.tsv"
LEIDEN_TAG = "leiden"
LINKAGES = ["average", "ward"]
MAPPERS = ["nearest_centroid", "wann", "tacco", "tangram", "dot"]

# Walltime classes (one array job each), from Experiment 1's slowest runs of these pairs
# (TinyFat, 8 cores) without the runs split below.
CLASSES = {
    "fast": "02:00:00",         # 0.5 h: nearest_centroid, wann, tacco
    "tangram": "02:00:00",      # 1.2 h
    "dot": "02:00:00",          # 1.6 h
}
# One job per K instead, then merged: every run of Experiment 1 from 1.8 h up.
SPLIT_BY_K = {
    (93, LEIDEN_TAG, "dot"),                                # 33.5 h
    (93, "annotation_cellType", "dot"),                     # 6.8 h
    (4, LEIDEN_TAG, "dot"),                                 # 14.3 h
    (4, "annotation_cell_ontology_class", "dot"),           # 9.2 h
    (4, "annotation_cell_ontology_relabel", "dot"),         # 4.4 h
    (2, LEIDEN_TAG, "dot"),                                 # 4.1 h
    (3, LEIDEN_TAG, "dot"),                                 # 3.8 h
    (94, "annotation_cellType", "dot"),                     # 3.7 h
    (94, LEIDEN_TAG, "dot"),                                # 2.6 h
    (93, LEIDEN_TAG, "tangram"),                            # 9.3 h
    (93, "annotation_cellType", "tangram"),                 # 3.1 h
    (4, LEIDEN_TAG, "tangram"),                             # 4.1 h
    (4, "annotation_cell_ontology_class", "tangram"),       # 2.9 h
    (4, "annotation_cell_ontology_relabel", "tangram"),     # 1.8 h
    (2, LEIDEN_TAG, "tangram"),                             # 3.4 h
}


def start_tag(ann_key: str) -> str:
    """The start folder name, as spatialaim.io.start_folder_name builds it."""
    return LEIDEN_TAG if not ann_key else "annotation_" + re.sub(r"[^A-Za-z0-9._-]+", "_", ann_key)


def data_dir() -> Path:
    return Path(os.environ.get("DEMO_DATA", HERE.parents[1] / "Data" / "04_Demo_sample" / "datasets"))


def runs() -> list[dict]:
    with open(RUNS, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def split_by_k(run: dict) -> bool:
    return (int(run["pair_id"]), run["start_tag"], run["mapper"]) in SPLIT_BY_K


def run_class(run: dict) -> str | None:
    """The walltime class of a sweep job; None for a run split by K."""
    if split_by_k(run):
        return None
    return run["mapper"] if run["mapper"] in CLASSES else "fast"


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


def write_runs(data: Path) -> None:
    """Rewrite runs.tsv from pairs.csv and the two index.csv files."""
    with open(data / "scRNA" / "index.csv", newline="", encoding="utf-8") as fh:
        keys = {row["Name"]: [k for k in ((row.get(f"CellTypeKey{i}") or "").strip() for i in range(3))
                              if k and k.lower() != "nan"]
                for row in csv.DictReader(fh)}
    with open(data / "ST" / "index.csv", newline="", encoding="utf-8") as fh:
        masks = {row["Name"]: (row.get("Segmentation") or "").strip() for row in csv.DictReader(fh)}
    with open(data / "pairs.csv", newline="", encoding="utf-8") as fh:
        pairs = sorted(csv.DictReader(fh), key=lambda p: int(p["PairID"]))
    rows = []
    for p in pairs:
        pair_dir = f"{int(p['PairID']):03d}_{p['scName']}__{p['stName']}"
        for ann_key in [""] + keys[p["scName"]]:
            for linkage in LINKAGES:
                for mapper in MAPPERS:
                    rows.append({"run": len(rows), "pair_id": int(p["PairID"]), "pair_dir": pair_dir,
                                 "sc": p["scName"], "st": p["stName"],
                                 "start_tag": start_tag(ann_key), "ann_key": ann_key,
                                 "linkage": linkage, "mapper": mapper,
                                 "segmentation": masks.get(p["stName"], "")})
    with open(RUNS, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    split = [r for r in rows if split_by_k(r)]
    print(f"{len(rows)} runs from {len(pairs)} pairs -> {RUNS} ({len(split)} split by K)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("runs", help="rewrite runs.tsv from $DEMO_DATA")
    sub.add_parser("pairs", help="print the scaffold jobs' --array list")
    for name in ("array", "time"):
        sub.add_parser(name).add_argument("cls", choices=list(CLASSES))
    args = ap.parse_args()
    if args.cmd == "runs":
        write_runs(data_dir())
    elif args.cmd == "pairs":
        print(array_spec({r["pair_id"] for r in runs()}))
    elif args.cmd == "time":
        print(CLASSES[args.cls])
    else:
        print(array_spec(r["run"] for r in runs() if run_class(r) == args.cls))


if __name__ == "__main__":
    main()
