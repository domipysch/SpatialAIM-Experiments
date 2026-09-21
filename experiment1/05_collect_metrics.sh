#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Experiment 1 - collect the metric files of every run root into one archive.
#
#   bash 05_collect_metrics.sh                 # collect + tar
#   bash 05_collect_metrics.sh --dry-run       # report only, stage nothing
#
# Run this unchanged on BOTH servers. The two archives hold disjoint mappers,
# share the same relative layout, and are merged by extracting both into the
# same folder.
#
# Collected per run root (<start_tag>/<pair_dir>/<tree>/<mapper>/):
#     k_comparison.csv    the metric curve over K - the whole of the analysis
#     k_selection.json    the proposal (best.overall)
#     config.yaml         provenance: which K range this sweep actually ran
# Plus, host-scoped under _collect/<host>/, the run ledgers and this run's
# manifest. No h5ad, no k_<kkk>/ folder, no log - the archive stays ~50 MB.
#
# Only the 92 WHOLE pairs are collected. Datasets 11 and 12 exist on disk both
# tiled (11_1_..., 12_1_...) and whole (11_human-breast-cancer,
# 12_human-mtg); the tiled run roots are skipped. The pair list comes from
# PAIRS_CSV when that file is present, otherwise from the name rule, and the
# script prints which of the two it used and how many pairs it kept - that
# number must be 92.
#
# Read-only on the results tree: everything is copied into STAGE first.
#
# Env vars:
#   OUT_ROOT   results root  (default /data/bionets/zi69hebi/Alignment/02_Results/experiment1)
#   PAIRS_CSV  whole-pair list
#              (default <OUT_ROOT>/../../01_Datasets/pairs_only_raw_counts_11_12_whole.csv)
#   STAGE      staging dir   (default $PWD/exp1_collect)
#   HOST_TAG   name for this host in _collect/ (default: hostname -s)
# ---------------------------------------------------------------------------
set -uo pipefail

OUT_ROOT="${OUT_ROOT:-/data/bionets/zi69hebi/Alignment/02_Results/experiment1}"
DATA_ROOT="${DATA_ROOT:-$(dirname "$(dirname "$OUT_ROOT")")/01_Datasets}"
PAIRS_CSV="${PAIRS_CSV:-$DATA_ROOT/pairs_only_raw_counts_11_12_whole.csv}"
STAGE="${STAGE:-$PWD/exp1_collect}"
HOST_TAG="${HOST_TAG:-$(hostname -s 2>/dev/null || echo host)}"
DRY_RUN=0
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=1

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

[[ -d "$OUT_ROOT" ]] || die "no such output root: $OUT_ROOT"
command -v python3 >/dev/null 2>&1 || die "python3 is not on PATH"
command -v tar     >/dev/null 2>&1 || die "tar is not on PATH"

ARCHIVE="$PWD/exp1_metrics_${HOST_TAG}_$(date +%Y%m%d).tar.gz"

printf '\n=== Experiment 1 metric collection\n'
printf '  results  : %s\n' "$OUT_ROOT"
printf '  pairs    : %s\n' "$PAIRS_CSV"
printf '  stage    : %s/experiment1\n' "$STAGE"
printf '  host tag : %s\n\n' "$HOST_TAG"

OUT_ROOT="$OUT_ROOT" PAIRS_CSV="$PAIRS_CSV" STAGE="$STAGE" \
HOST_TAG="$HOST_TAG" DRY_RUN="$DRY_RUN" python3 <<'PY' || die "collection failed"
import csv, os, re, shutil, sys
from collections import Counter
from pathlib import Path

out_root  = Path(os.environ["OUT_ROOT"])
pairs_csv = Path(os.environ["PAIRS_CSV"])
stage     = Path(os.environ["STAGE"]) / "experiment1"
host_tag  = os.environ["HOST_TAG"]
dry_run   = os.environ["DRY_RUN"] == "1"

WANTED = ("k_comparison.csv", "k_selection.json", "config.yaml")

# --- which pair directories count -------------------------------------------
# A pair dir is "<PairID:03d>_<scName>__<stName>". Preferred source of truth is
# the whole-pair CSV; the fallback rule drops only the tiled slices of datasets
# 11 and 12 (whose whole slices carry no slice index), and nothing else - the
# multi-slice datasets 05/07/09/10 all have indexed names and are all kept.
TILED = re.compile(r"__(11|12)_\d+_")

allowed, source = None, "name rule (11_<n>_/12_<n>_ dropped)"
if pairs_csv.exists():
    allowed = set()
    with open(pairs_csv, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            allowed.add(f"{int(row['PairID']):03d}_{row['scName'].strip()}"
                        f"__{row['stName'].strip()}")
    source = f"{pairs_csv.name} ({len(allowed)} pairs)"
else:
    print(f"  WARNING: {pairs_csv} not found - falling back to the name rule")

def keep(pair_dir: str) -> bool:
    return pair_dir in allowed if allowed is not None else not TILED.search(pair_dir)

# --- walk the run roots ------------------------------------------------------
n_roots = n_files = n_skipped_pairs = 0
per_mapper, per_stage_mode, missing = Counter(), Counter(), Counter()
skipped_pairs, kept_pairs = set(), set()
copied: list[str] = []

for table in sorted(out_root.rglob("k_comparison.csv")):
    run_root = table.parent
    try:
        parts = run_root.relative_to(out_root).parts
    except ValueError:
        continue
    if len(parts) != 4:                      # not <tag>/<pair>/<tree>/<mapper>
        continue
    start_tag, pair_dir, tree, mapper = parts
    if not keep(pair_dir):
        skipped_pairs.add(pair_dir)
        n_skipped_pairs += 1
        continue

    kept_pairs.add(pair_dir)
    n_roots += 1
    per_mapper[mapper] += 1
    per_stage_mode[("leiden" if start_tag.startswith("leiden") else "annotation", mapper)] += 1

    for name in WANTED:
        src = run_root / name
        if not src.is_file():
            missing[name] += 1
            continue
        n_files += 1
        rel = Path(*parts) / name
        copied.append(str(rel))
        if not dry_run:
            dst = stage / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)

# --- the host-scoped extras --------------------------------------------------
extras = 0
if not dry_run and n_roots:
    coll = stage / "_collect" / host_tag
    coll.mkdir(parents=True, exist_ok=True)
    for log_root in sorted(out_root.glob("_logs*")):
        for name in ("ledger.csv", "k_selection_index.csv"):
            src = log_root / name
            if src.is_file():
                shutil.copy2(src, coll / f"{log_root.name}__{name}")
                extras += 1
    (coll / "filelist.txt").write_text("\n".join(sorted(copied)) + "\n", encoding="utf-8")
    (coll / "manifest.txt").write_text(
        f"host      : {host_tag}\n"
        f"out_root  : {out_root}\n"
        f"pair list : {source}\n"
        f"run roots : {n_roots}\n"
        f"files     : {n_files}\n"
        f"mappers   : {', '.join(sorted(per_mapper))}\n"
        f"pairs     : {len(kept_pairs)}\n",
        encoding="utf-8",
    )

# --- report ------------------------------------------------------------------
print(f"  pair list : {source}")
print(f"  pairs kept: {len(kept_pairs)}   (expected 92)")
print(f"  run roots : {n_roots}   ({n_skipped_pairs} skipped in "
      f"{len(skipped_pairs)} tiled pair dir(s))")
print(f"  files     : {n_files}" + ("  [dry run - nothing copied]" if dry_run else ""))
if extras:
    print(f"  extras    : {extras} ledger/index file(s) under _collect/{host_tag}/")

print("\n  per mapper (expected 526 = 184 leiden + 342 annotation, when complete):")
for mapper in sorted(per_mapper):
    le = per_stage_mode[("leiden", mapper)]
    an = per_stage_mode[("annotation", mapper)]
    flag = "" if (le, an) == (184, 342) else "   <-- incomplete"
    print(f"    {mapper:<16} {per_mapper[mapper]:>4}   leiden {le:>3}/184   "
          f"annotation {an:>3}/342{flag}")

if missing:
    print("\n  missing files (run roots lacking them):")
    for name, count in sorted(missing.items()):
        print(f"    {name:<20} {count}")

if not n_roots:
    print("\n  nothing collected - check OUT_ROOT")
    sys.exit(1)
PY

if (( DRY_RUN )); then
  printf '\n--dry-run: nothing was staged or archived.\n\n'
  exit 0
fi

printf '\n  archiving -> %s\n' "$ARCHIVE"
tar czf "$ARCHIVE" -C "$STAGE" experiment1 || die "tar failed"
printf '  done: %s\n\n' "$(du -h "$ARCHIVE" | cut -f1)"
printf 'Copy it to the workstation with:\n  scp %s:%s .\n\n' \
  "$(hostname -s 2>/dev/null || echo '<server>')" "$ARCHIVE"