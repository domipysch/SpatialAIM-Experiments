"""
Step 1 of figure 4's kidney example: compute the two mappings the figure shows.

Experiment 1 compares, per pair, mapper and tree, SpatialAIM's mapping with the
mapping of the author's annotation. This script reads which two mappings those
are from the Experiment 1 summary and recomputes exactly them -- the per-spot
results of Experiment 1 stay on the cluster, only its summary is pulled:

    annotation   the annotation start at K = T (the identity cut, i.e. every spot
                 mapped straight onto the author's types); where a reference has
                 two annotations, the one exp1.paired compares against
    spatialaim   the start exp1.best_start picks (Leiden or an annotation), at
                 its K = best.overall

Each is one ``spatialaim run`` with --k_min = --k_max = that K and only that
side's start (an annotation start with --no_leiden_start), so no sweep is
repeated. The two pair outputs go to

    Data/02_Results/biology_example/<pair_dir>/<mapper>_<tree>/{annotation,spatialaim}/

each holding the run root ``<start>/<mapper>/`` (``leiden`` or ``annotation_<key>``).

next to a ``selection.json`` recording what was chosen and why. Afterwards the
metrics of both runs are compared with the summary's; they should agree up to the
small run-to-run noise of the Leiden over-clustering.

    python sections/04_biology/04_kidney_example/01_run_mappings.py [--pair 4] [--ann_key KEY] [--force]
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
from spatialaim.io import start_folder_name

from spatialaim_experiments import exp1, paths
from spatialaim_experiments.style import MAPPERS

logger = logging.getLogger(__name__)

# Experiment 1 over-clustered every reference at this Leiden resolution.
LEIDEN_RESOLUTION = 3.0

# The metrics compared with the summary after the runs.
METRICS = ["spatial_assort_celltype", "spatial_assort_global",
           "modularity_st_expression", "modularity_shared"]

OUT = paths.RESULTS / "biology_example"


# ---------------------------------------------------------------------------
# which two mappings
# ---------------------------------------------------------------------------

def select(pair: int, mapper: str, tree: str, ann_key: str | None = None) -> dict:
    """The two mappings Experiment 1 compares for one pair, mapper and tree;
    ``ann_key`` swaps the annotation side for another of the reference's annotations
    (at its own K = T), e.g. the kidney's cell_ontology_relabel.

    Output
        {"pair_id", "pair_dir", "reference", "st", "mapper", "tree",
         "annotation": {"ann_key", "k", <metrics>},
         "spatialaim": {"ann_key" (None = Leiden start), "k", <metrics>}}
    """
    df = exp1.best_start(exp1.load())
    rows = df[(df["pair_id"] == pair) & (df["mapper"] == mapper) & (df["tree"] == tree)]
    if rows.empty:
        raise SystemExit(f"no Experiment 1 runs for pair {pair}, {mapper}, {tree}")

    # exactly the pairing the figure's scatters use, unless ann_key names another
    # of the reference's annotations
    paired = exp1.paired(rows, METRICS[0]).iloc[0]
    key = ann_key or paired["ann_key_base"]
    spatialaim = rows[rows["role"] == "spatialaim"].iloc[0]
    baseline = rows[(rows["role"] == "baseline") & (rows["ann_key"] == key)]
    if baseline.empty:
        raise SystemExit(f"no Experiment 1 baseline with annotation {key!r} for pair {pair}")
    baseline = baseline.iloc[0]

    pairs = pd.read_csv(paths.PAIRS_CSV)
    st = pairs.loc[pairs["PairID"] == pair, "stName"].iloc[0]

    def side(row) -> dict:
        key = row["ann_key"]
        return {"ann_key": None if pd.isna(key) else str(key), "k": int(row["k"]),
                **{m: float(row[m]) for m in METRICS}}

    return {"pair_id": pair, "pair_dir": spatialaim["pair_dir"],
            "reference": spatialaim["reference"], "st": st,
            "mapper": mapper, "tree": tree,
            "annotation": side(baseline), "spatialaim": side(spatialaim)}


# ---------------------------------------------------------------------------
# running them
# ---------------------------------------------------------------------------

def conda_env() -> dict:
    """The environment for the ``spatialaim run`` subprocess: this Python's conda
    env first on PATH. Without its DLL folders the squidpy import in the
    analysis step crashes a fresh subprocess on Windows."""
    prefix = Path(sys.prefix)
    dirs = [prefix, prefix / "Library" / "mingw-w64" / "bin", prefix / "Library" / "usr" / "bin",
            prefix / "Library" / "bin", prefix / "Scripts", prefix / "bin"]
    env = dict(os.environ)
    env["PATH"] = os.pathsep.join([*map(str, dirs), env.get("PATH", "")])
    return env


def pair_output(selection: dict, side: str) -> Path:
    """Where one side's ``spatialaim run`` writes (its ``--output_dir``)."""
    return (OUT / selection["pair_dir"] / f"{selection['mapper']}_{selection['tree']}"
            / side)


def run_root(selection: dict, side: str) -> Path:
    """One side's run root, ``<pair output>/<start>/<mapper>``; the start folder
    above it holds ``start_clustering.h5ad``."""
    key = selection[side]["ann_key"]
    start = "leiden" if key is None else start_folder_name(key)
    return pair_output(selection, side) / start / selection["mapper"]


def run(selection: dict, side: str, force: bool) -> Path:
    """One ``spatialaim run`` at the single K this side needs; skipped when that
    K's mapping exists already (unless ``force``)."""
    spec = selection[side]
    out = pair_output(selection, side)
    k_dir = run_root(selection, side) / f"k_{spec['k']:03d}"
    if k_dir.exists() and not force:
        logger.info("%-10s exists, skipped: %s", side, k_dir)
        return out

    cmd = [sys.executable, "-m", "spatialaim", "run",
           "--scdata", str(paths.DATASETS / "scRNA" / f"{selection['reference']}.h5ad"),
           "--stdata", str(paths.DATASETS / "ST" / f"{selection['st']}.h5ad"),
           "--output_dir", str(out),
           "--mapping", selection["mapper"],
           "--linkage_method", selection["tree"],
           "--leiden_resolution", str(LEIDEN_RESOLUTION),
           "--k_min", str(spec["k"]), "--k_max", str(spec["k"])]
    if spec["ann_key"] is not None:
        cmd += ["--annotation", spec["ann_key"], "--no_leiden_start"]
    logger.info("%-10s running at K = %d ...", side, spec["k"])
    subprocess.run(cmd, check=True, env=conda_env())
    return out


def compare(selection: dict, side: str) -> None:
    """Log this side's metrics next to Experiment 1's."""
    spec = selection[side]
    table = pd.read_csv(run_root(selection, side) / "k_comparison.csv")
    row = table[table["k"] == spec["k"]].iloc[0]
    for metric in METRICS:
        logger.info("%-10s %-26s here %10.3f   Experiment 1 %10.3f",
                    side, metric, row[metric], spec[metric])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pair", type=int, default=4, help="PairID (default: 4, the kidney)")
    ap.add_argument("--mapper", choices=MAPPERS, default="nearest_centroid")
    ap.add_argument("--tree", choices=["average", "ward"], default="average")
    ap.add_argument("--ann_key", default=None,
                    help="annotation side: this obs column instead of the one exp1.paired picks")
    ap.add_argument("--force", action="store_true", help="recompute existing mappings")
    args = ap.parse_args()

    selection = select(args.pair, args.mapper, args.tree, args.ann_key)
    folder = pair_output(selection, "annotation").parent
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "selection.json").write_text(json.dumps(selection, indent=2), encoding="utf-8")
    logger.info("pair %d (%s): annotation %s at K = %d, SpatialAIM %s start at K = %d",
                selection["pair_id"], selection["reference"],
                selection["annotation"]["ann_key"], selection["annotation"]["k"],
                selection["spatialaim"]["ann_key"] or "Leiden", selection["spatialaim"]["k"])

    for side in ("annotation", "spatialaim"):
        run(selection, side, args.force)
        compare(selection, side)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
