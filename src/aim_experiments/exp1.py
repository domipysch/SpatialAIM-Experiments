"""Experiment 1 -> one tidy table, shared by the reconstruction and biology
figures (sections 03 and 04).

Reads only the three small files collected from the servers
(k_comparison.csv, k_selection.json, config.yaml) at each of the 2630 run roots
under ``Data/02_Results/experiment1``, and reduces each to the two rows the
figures need:

    baseline   annotation mode, K = T (the identity cut: every spot maps
               straight onto the author's annotation)
    aim        Leiden mode, K = best.overall from k_selection.json

Both are read from the same column of the same metric table, so the comparison
is apples to apples.

Two traps this handles:

* Windows merged `ann_cellType` and `ann_celltype` into one directory, so the
  start tag in the path is not the annotation key. The true key is read from
  each run root's config.yaml (`start_from_annotation`).
* Only the 92 whole pairs count; the tiled 11_<n>/12_<n> run roots are not in
  the collected tree, but the pair list is applied anyway so this cannot drift.

Aggregation is always per reference, never pooled over pairs: 71 of the 92
pairs are retina islands.

Usage
-----
    python -m aim_experiments.exp1            # build + cache + summarise
    python -m aim_experiments.exp1 --force    # ignore the cache
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import logging
import sys
from pathlib import Path

import pandas as pd
import yaml

from . import paths
from .style import MAPPERS, TREES  # noqa: F401  (re-exported: callers use exp1.MAPPERS)

logger = logging.getLogger(__name__)


def _kselection():
    """``aim.metrics.kselection`` from the installed ``spatial-aim``, loaded
    straight from its file.

    Importing it by name would first execute ``aim/metrics/__init__.py``, which
    pulls in the scanpy/anndata stack for metrics this module never touches. The
    module itself needs only numpy and pandas. Its *location* still comes from
    the installed package -- ``aim/__init__.py`` is import-light by design -- so
    this is the same scoring code that wrote every ``k_selection.json`` on disk,
    without paying for the rest of the package.
    """
    import aim

    path = Path(aim.__file__).parent / "metrics" / "kselection.py"
    spec = importlib.util.spec_from_file_location("aim_kselection", path)
    module = importlib.util.module_from_spec(spec)
    # @dataclass resolves its own module through sys.modules, so register first
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


KS = _kselection()

CACHE = paths.outputs("experiment1") / "runs.csv"

# The two reconstruction curves section 03 is about: hard mapping, raw counts.
METRICS = ["cossim_hard_raw_spot", "cossim_hard_raw_gene"]
NULLS = [f"{m}_null" for m in METRICS]

# The other two criteria's curves, carried along for section 04 (biology and
# topology). They have no null column: the spatial pair are z-scores against the
# label shuffle already, and modularity is defined against a degree-preserving
# random partition, so chance is 0 for all four by construction.
BIO_METRICS = [
    "nhood_mean_self_zscore",
    "local_purity_zscore",
    "modularity_st_expression",
    "modularity_shared",
]


# ---------------------------------------------------------------------------
# the corpus
# ---------------------------------------------------------------------------

def pairs() -> pd.DataFrame:
    """The 92 whole pairs, with the reference's annotation keys attached."""
    df = pd.read_csv(paths.PAIRS_CSV)
    sc = pd.read_csv(paths.SC_INDEX).set_index("Name")
    keys = {}
    for name, row in sc.iterrows():
        found = [str(row.get(f"CellTypeKey{i}", "")).strip() for i in range(3)]
        keys[name] = [k for k in found if k and k.lower() != "nan"]
    df["ann_keys"] = df["scName"].map(keys)
    df["pair_dir"] = [
        f"{int(r.PairID):03d}_{r.scName}__{r.stName}" for r in df.itertuples()
    ]
    # key0 is the author's primary annotation and the baseline used in the
    # AIM-vs-annotation panel; key1, where it exists, is the second granularity.
    df["key0"] = [k[0] if k else None for k in df["ann_keys"]]
    df["key1"] = [k[1] if len(k) > 1 else None for k in df["ann_keys"]]
    return df


# ---------------------------------------------------------------------------
# one row per run root
# ---------------------------------------------------------------------------

def _read_run_root(run_root: Path, out_root: Path) -> dict | None:
    """The one K-level that matters for this run root, plus its context."""
    start_tag, pair_dir, tree, mapper = run_root.relative_to(out_root).parts
    table = pd.read_csv(run_root / "k_comparison.csv")
    if "k" not in table.columns or table.empty:
        return None
    table = table.set_index(table["k"].astype(int))

    cfg = yaml.safe_load((run_root / "config.yaml").read_text(encoding="utf-8")) or {}
    ann_key = cfg.get("start_from_annotation")     # None in Leiden mode

    n_start = int(table.index.max())               # T in annotation mode, L in Leiden
    if ann_key is None:
        selection = json.loads((run_root / "k_selection.json").read_text(encoding="utf-8"))
        best = (selection.get("best") or {}).get("overall")
        if best is None:
            return None
        k, role = int(best["k"]), "aim"
    else:
        # the identity cut: the author's annotation, mapped as it stands
        k, role = n_start, "baseline"

    if k not in table.index:
        return None
    row = table.loc[k]

    # The overall score at the level that is actually used -- AIM's proposal, or
    # the annotation's identity cut. Recorded for *both* roles because it is what
    # decides which of a reference's annotations becomes the baseline: AIM's K is
    # chosen by the overall score, so the annotation must be too, and not by
    # whichever measure a given panel happens to plot.
    scores = KS.score_table(table.reset_index(drop=True))
    overall = pd.Series(KS.overall(scores), index=scores["k"].astype(int))
    score = float(overall.get(k, float("nan")))

    record = {
        "role": role,
        "start_tag": start_tag,
        "ann_key": ann_key or "",
        "pair_dir": pair_dir,
        "tree": tree,
        "mapper": mapper,
        "k": k,
        "n_start": n_start,
        "score_overall": score,
        # the K range actually swept, so a truncated sweep cannot pass unnoticed
        "k_min_swept": int(table.index.min()),
        "n_k": int(len(table)),
    }
    for column in METRICS + NULLS + BIO_METRICS:
        record[column] = float(row[column]) if column in table.columns else float("nan")
    return record


def build(out_root: Path | None = None) -> pd.DataFrame:
    """Walk every run root once. ~2630 small CSV reads, about a minute."""
    out_root = paths.EXPERIMENT1 if out_root is None else Path(out_root)
    keep = set(pairs()["pair_dir"])
    rows, skipped = [], 0
    for table_path in sorted(out_root.rglob("k_comparison.csv")):
        run_root = table_path.parent
        parts = run_root.relative_to(out_root).parts
        if len(parts) != 4 or parts[1] not in keep:
            skipped += 1
            continue
        record = _read_run_root(run_root, out_root)
        if record is None:
            skipped += 1
            continue
        rows.append(record)
    logger.info("%d run roots read, %d skipped", len(rows), skipped)

    df = pd.DataFrame(rows)
    meta = pairs().set_index("pair_dir")
    df["pair_id"] = df["pair_dir"].map(meta["PairID"]).astype(int)
    df["reference"] = df["pair_dir"].map(meta["scName"])
    df["slice"] = df["pair_dir"].map(meta["stName"])
    df["key0"] = df["pair_dir"].map(meta["key0"])
    df["key1"] = df["pair_dir"].map(meta["key1"])
    # which granularity this baseline is: the author's primary, or the second one
    df["granularity"] = "-"
    df.loc[df["ann_key"] == df["key0"], "granularity"] = "key0"
    df.loc[df["ann_key"] == df["key1"], "granularity"] = "key1"
    return df.sort_values(["pair_id", "start_tag", "tree", "mapper"]).reset_index(drop=True)


def load(force: bool = False) -> pd.DataFrame:
    """The cached table, built on first use. Delete ``CACHE`` or pass ``force``
    after new runs have been collected."""
    if CACHE.exists() and not force:
        return pd.read_csv(CACHE)
    df = build()
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(CACHE, index=False)
    logger.info("cached -> %s", CACHE)
    return df


# ---------------------------------------------------------------------------
# the two comparisons the figures make
# ---------------------------------------------------------------------------

def paired(df: pd.DataFrame, metric: str = "cossim_hard_raw_spot",
           baseline: str = "best") -> pd.DataFrame:
    """AIM (Leiden, proposed K) against the annotation baseline (K = T), joined
    per (pair, tree, mapper). One row per comparison; `delta` is AIM - baseline.

    ``baseline`` picks which annotation AIM is measured against where the
    reference carries more than one:

    ``best``  the one whose identity cut scores highest **overall** -- the
              harmonic mean across reconstruction, spatial organisation and
              coherence, exactly the score by which AIM picks its own K. The
              default, and the only symmetric choice: selecting the annotation on
              the measure a panel happens to plot would hand the baseline a
              per-panel advantage AIM is never given.
    ``key0``  the annotation the dataset index lists first.

    Which of two annotations is called the primary one in the dataset index is
    arbitrary, and it matters: on 01_mouse-ssp, AIM beats ``cellType`` (3 types)
    by +0.26 on per-spot cosine and loses to ``cellTypeMinor`` (22 types) by
    -0.05. Comparing against the better of the author's own options removes that
    lever.

    Both sides always use the same mapper, so the difference is the starting
    partition, not the aligner.
    """
    aim = df[df["role"] == "aim"]
    base = df[df["role"] == "baseline"]
    join = ["pair_dir", "tree", "mapper"]
    if baseline == "key0":
        base = base[base["granularity"] == "key0"]
    elif baseline == "best":
        base = base.loc[base.groupby(join)["score_overall"].idxmax()]
    else:
        raise ValueError(f"baseline must be 'best' or 'key0', got {baseline!r}")

    out = aim.merge(
        base[join + [metric, "k", "ann_key", "granularity"]],
        on=join, suffixes=("", "_base"), validate="one_to_one",
    )
    out["delta"] = out[metric] - out[f"{metric}_base"]
    out["k_ratio"] = out["k"] / out["k_base"]
    return out


def granularity(df: pd.DataFrame, metric: str = "cossim_hard_raw_spot") -> pd.DataFrame:
    """The two author granularities against each other, per (pair, tree, mapper).

    Only the five references that carry a second annotation key contribute.
    """
    base = df[df["role"] == "baseline"]
    join = ["pair_dir", "tree", "mapper"]
    coarse = base[base["granularity"] == "key0"]
    fine = base[base["granularity"] == "key1"]
    out = coarse.merge(
        fine[join + [metric, "k", "ann_key"]],
        on=join, suffixes=("_key0", "_key1"), validate="one_to_one",
    )
    out["delta"] = out[f"{metric}_key1"] - out[f"{metric}_key0"]
    return out


def per_reference(frame: pd.DataFrame, value: str = "delta") -> pd.DataFrame:
    """Collapse pairs to their reference -- the only aggregation this corpus
    permits, since 71 of 92 pairs are retina islands."""
    return (
        frame.groupby(["reference", "mapper", "tree"], as_index=False)
        .agg(mean=(value, "mean"), sd=(value, "std"), n=(value, "size"))
    )


# ---------------------------------------------------------------------------

def summarise(df: pd.DataFrame) -> None:
    print(f"\nrun roots: {len(df)}   references: {df['reference'].nunique()}   "
          f"pairs: {df['pair_dir'].nunique()}")
    print("\nroles:")
    print(df.groupby(["role", "granularity"]).size().to_string())

    truncated = df[df["k_min_swept"] > 1]
    print(f"\nsweeps not reaching K=1: {len(truncated)}")

    for metric in METRICS:
        pr = paired(df, metric)
        agg = per_reference(pr)
        wins = (agg["mean"] > 0).sum()
        print(f"\n=== {metric}")
        print(f"  AIM vs annotation: {len(pr)} comparisons, "
              f"{(pr['delta'] > 0).mean():.1%} positive; "
              f"per (reference,mapper,tree) cells positive: {wins}/{len(agg)}")
        print("  mean delta per reference (averaged over mappers and trees):")
        print(pr.groupby("reference")["delta"].agg(["mean", "size"]).round(4).to_string())

        gr = granularity(df, metric)
        if len(gr):
            print(f"  granularity key0 -> key1: {len(gr)} comparisons")
            print(gr.groupby("reference")["delta"].agg(["mean", "size"]).round(4).to_string())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true", help="rebuild, ignore the cache")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    summarise(load(force=args.force))


if __name__ == "__main__":
    main()
