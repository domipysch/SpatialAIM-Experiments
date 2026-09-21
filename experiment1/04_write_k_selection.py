#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Experiment 1 - write the missing k_selection.json for every run root.
#
#   conda run -n spatialaim_env python 04_write_k_selection.py
#   conda run -n spatialaim_env python 04_write_k_selection.py --dry-run
#
# The sweeps were run with a version that wrote k_comparison.csv but no
# k_selection.json. Nothing is lost: the proposal is derived from that table
# alone (spatialaim.analysis.ksweep.select_k reads no other file), so this recomputes
# it for every run root without touching a single mapping.
#
# It also detects *truncated* tables. A sweep resumed with --k_max writes a
# k_comparison.csv covering only the resumed sub-range, while the K levels of
# the earlier attempt still sit on disk. Those roots are reported; --repair
# additionally rebuilds them with compare_k_runs, which only re-reads each
# level's analysis/data, so no alignment is re-run there either.
#
# Side output: one consolidated CSV of every proposal (--index), handy as a
# coverage check before pulling the metrics off the server.
#
# What it writes, and nothing else:
#   <run_root>/k_selection.json          created (this is the point)
#   <run_root>/k_comparison.csv          overwritten ONLY with --repair, ONLY at
#                                        truncated roots, ONLY over a superset of
#                                        the rows already there
#   <out_root>/_logs/k_selection_index.csv   created
# Every k_<kkk>/ folder, every analysis/data file, every h5ad and every log is
# read at most.
# ---------------------------------------------------------------------------

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path

import pandas as pd

from spatialaim.analysis.ksweep import compare_k_runs, select_k

DEFAULT_OUT_ROOT = "/data/bionets/zi69hebi/Alignment/02_Results/experiment1"

# What score_table() consumes. A table missing any of these yields an all-NaN
# criterion and therefore no proposal at all, so say so loudly rather than
# writing 2630 empty JSONs.
REQUIRED = (
    "k",
    "cossim_hard_raw_spot",
    "cossim_hard_raw_gene",
    "nhood_mean_self_zscore",
    "local_purity_zscore",
    "modularity_st_expression",
    "modularity_shared",
)


def completed_ks(run_root: Path) -> list[int]:
    """The K levels that actually finished their analysis, from the k_<kkk>/ dirs."""
    return sorted(
        int(p.name[2:])
        for p in run_root.glob("k_[0-9]*")
        if (p / "analysis" / "data" / "cossim_summary.csv").is_file()
    )


def parse_run_root(run_root: Path, out_root: Path) -> dict:
    """Split <start_tag>/<PairID>_<sc>__<st>/<tree>/<mapper> back into fields."""
    parts = run_root.relative_to(out_root).parts
    if len(parts) != 4:
        return {}
    start_tag, pair_dir, tree, mapper = parts
    pair_id, _, rest = pair_dir.partition("_")
    sc, _, st = rest.partition("__")
    return {
        "start_tag": start_tag,
        "pair_id": pair_id,
        "sc": sc,
        "st": st,
        "tree": tree,
        "mapper": mapper,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-root", default=DEFAULT_OUT_ROOT)
    ap.add_argument(
        "--index",
        default="",
        help="write a consolidated CSV of all proposals here "
        "(default: <out-root>/_logs/k_selection_index.csv)",
    )
    ap.add_argument(
        "--repair",
        action="store_true",
        help="also rebuild k_comparison.csv where the table is missing completed "
        "K levels (overwrites that file). Off by default, so a plain run only "
        "ever creates k_selection.json; truncated roots are reported either way.",
    )
    ap.add_argument("--dry-run", action="store_true", help="report, write nothing")
    args = ap.parse_args()

    # select_k / compare_k_runs log one INFO line per call; 2630 of those bury
    # the summary that matters.
    logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")

    out_root = Path(args.out_root)
    if not out_root.is_dir():
        print(f"ERROR: no such output root: {out_root}", file=sys.stderr)
        return 1

    tables = sorted(out_root.rglob("k_comparison.csv"))
    print(f"{len(tables)} run root(s) under {out_root}")
    if not tables:
        return 1

    index_path = Path(args.index) if args.index else out_root / "_logs" / "k_selection_index.csv"

    rows: list[dict] = []
    n_written = n_repaired = n_no_proposal = 0
    problems: list[str] = []

    for i, table_path in enumerate(tables, 1):
        run_root = table_path.parent
        rel = run_root.relative_to(out_root)
        try:
            df = pd.read_csv(table_path)
        except Exception as exc:
            problems.append(f"unreadable table  {rel}: {exc}")
            continue

        missing = [c for c in REQUIRED if c not in df.columns]
        if missing:
            problems.append(f"missing columns  {rel}: {', '.join(missing)}")
            continue

        # Truncated table? The k_<kkk>/ folders are the ground truth for what ran.
        on_disk = completed_ks(run_root)
        in_table = sorted(int(k) for k in df["k"])
        if on_disk and set(on_disk) - set(in_table):
            extra = len(set(on_disk) - set(in_table))
            if not args.repair:
                problems.append(
                    f"truncated table   {rel}: {extra} completed K level(s) not in "
                    f"the table (pass --repair to rebuild it)"
                )
            elif args.dry_run:
                problems.append(f"would rebuild     {rel}: +{extra} K level(s)")
                n_repaired += 1
            else:
                # Rebuild over the *union*, never over the on-disk levels alone: a
                # row the old table had must not vanish because its k_<kkk>/ folder
                # is gone. Levels without analysis data come back as NaN, exactly
                # as collect_k_metrics produced them the first time.
                compare_k_runs(run_root, sorted(set(on_disk) | set(in_table)))
                df = pd.read_csv(table_path)
                n_repaired += 1

        if args.dry_run:
            continue

        try:
            proposal = select_k(run_root, df)       # writes k_selection.json
        except Exception as exc:
            problems.append(f"select_k failed   {rel}: {exc}")
            continue
        n_written += 1

        best = proposal.get("best", {})
        overall = best.get("overall")
        if overall is None:
            n_no_proposal += 1
            problems.append(f"no finite score   {rel}")

        row = parse_run_root(run_root, out_root)
        row.update(
            {
                "k_min_swept": in_table[0] if in_table else "",
                "k_max_swept": in_table[-1] if in_table else "",
                "n_k": len(in_table),
                "k_overall": overall["k"] if overall else "",
                "score_overall": f"{overall['score']:.6f}" if overall else "",
            }
        )
        for key in ("recon", "spatial", "coherence"):
            hit = best.get(key)
            row[f"k_{key}"] = hit["k"] if hit else ""
        rows.append(row)

        if i % 250 == 0:
            print(f"  ... {i}/{len(tables)}")

    if rows and not args.dry_run:
        index_path.parent.mkdir(parents=True, exist_ok=True)
        with open(index_path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    print()
    print(f"  k_selection.json written : {n_written}")
    print(f"  k_comparison.csv rebuilt : {n_repaired}")
    print(f"  no finite proposal       : {n_no_proposal}")
    if rows and not args.dry_run:
        print(f"  index                    : {index_path}")
    if problems:
        print(f"\n  PROBLEMS ({len(problems)}):")
        for line in problems[:40]:
            print(f"    {line}")
        if len(problems) > 40:
            print(f"    ... and {len(problems) - 40} more")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
