#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Experiment 1 - the full prerequisite grid, run by WORKERS shell workers.
#
#   pairs in pairs_only_raw_counts.csv
#     x start clusters { Leiden over-clustering } u { every cell-type annotation
#                        the pair's scRNA reference has, per scRNA/index.csv }
#     x agglomeration tree { ward, average }
#     x spot->state mapper { nearest_centroid, wann, tangram, tacco, dot }
#
# One `aim run` (a full K sweep) per combination. Resumable: a combination whose
# k_comparison.csv already exists is skipped, so the script can be killed and
# restarted at any time.
#
# Scheduling - fully decided up front, no coordination while running:
#   * unit of work ("group") = one (pair, start mode); it carries both trees, so a
#     group's reference scaffold is built once and hardlinked ward -> average
#   * the group list is deduplicated and then split across the workers before
#     anything starts (greedy longest-processing-time-first on CellCount x
#     SpotCount, so the workers get comparable total load). Each worker gets its
#     own work file in _logs/work/ and only ever touches those groups - two workers
#     can therefore never land on the same run root.
#   * one stage per mapper, in MAPPINGS order, with a barrier between stages: the
#     cheap in-process mappers finish across the whole grid first, and every later
#     stage reads the scaffold cache the first stage wrote
#
# Output layout (a run root is named after the mapper alone, so the start mode and
# the tree method both have to be directory levels above it):
#
#   $OUT_ROOT/<start_tag>/<PairID>_<scName>__<stName>/<tree>/<mapper>/k_<kkk>/...
#            <start_tag> = leiden_res3.0 | ann_<obs_column>
#
# Usage:
#   WORKERS=8 bash 02_run_experiment1.sh        # 8 sweeps at a time
#   DRY_RUN=1 bash 02_run_experiment1.sh        # preflight: split + counts only
#   MAPPINGS="nearest_centroid wann" WORKERS=8 bash 02_run_experiment1.sh
#
# Env vars (all optional):
#   WORKERS        parallel workers                          (default 1)
#   DATA_ROOT      dataset root         (default /data/bionets/zi69hebi/Alignment/01_Datasets)
#   OUT_ROOT       results root         (default <DATA_ROOT>/../02_Results/experiment1)
#   PAIRS_CSV      pairs file           (default $DATA_ROOT/pairs_only_raw_counts.csv)
#   AIM_ENV        conda env for AIM    (default aim_env)
#   MAPPINGS       mapper list; also the stage order
#   TREES          "ward average"
#   START_MODES    "leiden annotation"  (drop one to run only the other)
#   LEIDEN_RES     --leiden_resolution  (default 3.0)
#   K_MIN/K_MAX/K_STEP  passed through when set (default: full sweep, step 1)
#   JOB_TIMEOUT    e.g. 12h - kill a single sweep that hangs (default: no limit)
#   CACHE_ROOT     matplotlib/numba/XDG caches, kept off $HOME
#                                       (default <DATA_ROOT>/../.cache)
#   CONFIG_ROOT    XDG_CONFIG_HOME      (default <DATA_ROOT>/../.config)
#   OMP_NUM_THREADS  per worker (default: nproc / WORKERS, at least 1)
#   ORDER          lpt (default) = weighted greedy split, biggest groups first;
#                  file = keep pairs.csv order and deal the groups round-robin
#   SHARD_COUNT/SHARD_INDEX   split across machines by PairID % SHARD_COUNT
#   FORCE=1        start even though another instance of this script is running
#   PRUNE_SCAFFOLD=1  delete reference_scaffold.h5ad after the last mapper stage of
#                     a group - saves a lot of disk; the GUI rebuilds it on demand
#   LINK_SCAFFOLD=0   do not hardlink the ward scaffold cache into the average cell
#   DRY_RUN=1         enumerate + check annotations, run nothing
# ---------------------------------------------------------------------------
set -uo pipefail

DATA_ROOT="${DATA_ROOT:-/data/bionets/zi69hebi/Alignment/01_Datasets}"
OUT_ROOT="${OUT_ROOT:-$(dirname "$DATA_ROOT")/02_Results/experiment1}"
PAIRS_CSV="${PAIRS_CSV:-$DATA_ROOT/pairs_only_raw_counts.csv}"
SC_DIR="$DATA_ROOT/scRNA"
ST_DIR="$DATA_ROOT/ST"
SC_INDEX="$SC_DIR/index.csv"
ST_INDEX="$ST_DIR/index.csv"

AIM_ENV="${AIM_ENV:-aim_env}"
WORKERS="${WORKERS:-1}"
MAPPINGS="${MAPPINGS:-nearest_centroid wann tangram tacco dot}"
TREES="${TREES:-ward average}"
START_MODES="${START_MODES:-leiden annotation}"
LEIDEN_RES="${LEIDEN_RES:-3.0}"
K_MIN="${K_MIN:-}"
K_MAX="${K_MAX:-}"
K_STEP="${K_STEP:-}"
JOB_TIMEOUT="${JOB_TIMEOUT:-}"
SHARD_COUNT="${SHARD_COUNT:-1}"
SHARD_INDEX="${SHARD_INDEX:-0}"
ORDER="${ORDER:-lpt}"
FORCE="${FORCE:-0}"
PRUNE_SCAFFOLD="${PRUNE_SCAFFOLD:-0}"
LINK_SCAFFOLD="${LINK_SCAFFOLD:-1}"
DRY_RUN="${DRY_RUN:-0}"

LOG_ROOT="${LOG_ROOT:-$OUT_ROOT/_logs}"
LEDGER="$LOG_ROOT/ledger.csv"
WORK_DIR="$LOG_ROOT/work"                     # one work file per worker
GROUPS_FILE="$WORK_DIR/groups_assigned.tsv"   # the whole split, worker id in col 1
NO_ANN="-"          # sentinel: this group uses the Leiden over-clustering

log() { printf '%s | %s\n' "$(date '+%F %T')" "$*"; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

# --- environment ------------------------------------------------------------
command -v conda >/dev/null 2>&1 || die "conda is not on PATH"
CONDA_BASE="$(conda info --base)" || die "conda info --base failed"
# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"
export CONDA_EXE="${CONDA_EXE:-$CONDA_BASE/bin/conda}"   # AIM needs to find conda

# --- keep every cache off $HOME ---------------------------------------------
# The 2026-08-12 run lost 582 of 634 sweeps to this: the home directory went
# unwritable at 22:54 and from then on every single `conda run -n tangram_env ...`
# died before Python started ("An unexpected error has occurred") because conda
# writes into ~/.conda; matplotlib hit the same wall one line earlier and only
# warned. Nothing here belongs on the home filesystem - conda's package cache and
# envs live under $CONDA_BASE, and the tool caches can live next to the data.
CACHE_ROOT="${CACHE_ROOT:-$(dirname "$DATA_ROOT")/.cache}"
CONFIG_ROOT="${CONFIG_ROOT:-$(dirname "$DATA_ROOT")/.config}"
export MPLCONFIGDIR="$CACHE_ROOT/mpl"
export XDG_CACHE_HOME="$CACHE_ROOT"
export XDG_CONFIG_HOME="$CONFIG_ROOT"
export CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-$CONDA_BASE/pkgs}"
export CONDA_ENVS_DIRS="${CONDA_ENVS_DIRS:-$CONDA_BASE/envs}"
export NUMBA_CACHE_DIR="$CACHE_ROOT/numba"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME" "$XDG_CONFIG_HOME" "$NUMBA_CACHE_DIR" \
  || die "cannot create the cache dirs under $CACHE_ROOT"
for d in "$MPLCONFIGDIR" "$XDG_CACHE_HOME" "$XDG_CONFIG_HOME" "$NUMBA_CACHE_DIR" \
         "$CONDA_PKGS_DIRS" "$CONDA_ENVS_DIRS"; do
  [[ -w "$d" ]] || die "not writable: $d (fix this before starting - it is what killed the last run)"
done
# $HOME is no longer on the hot path, but conda still reads ~/.condarc and a full
# filesystem there breaks other things; say so rather than silently continuing.
[[ -w "${HOME:-/nonexistent}" ]] || printf 'WARNING: %s is not writable - caches are redirected to %s, continuing\n' \
  "${HOME:-\$HOME}" "$CACHE_ROOT" >&2

conda activate "$AIM_ENV" || die "cannot activate $AIM_ENV"
command -v aim >/dev/null 2>&1 || die "'aim' not found in $AIM_ENV - run 01_setup_envs.sh"

[[ "$WORKERS" =~ ^[0-9]+$ ]] && (( WORKERS >= 1 )) || die "WORKERS must be a positive integer"

# Split the node's cores across the workers; a single sweep is one process, so
# without this every worker's BLAS/numba would try to use the whole machine.
NPROC="$( { nproc || getconf _NPROCESSORS_ONLN; } 2>/dev/null || echo 8 )"
if [[ -z "${OMP_NUM_THREADS:-}" ]]; then
  OMP_NUM_THREADS=$(( NPROC / WORKERS ))
  (( OMP_NUM_THREADS >= 1 )) || OMP_NUM_THREADS=1
fi
export OMP_NUM_THREADS
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
export NUMBA_NUM_THREADS="$OMP_NUM_THREADS"

# Safety net, not the fix: the aligner wrappers now build their h5ad indices with an
# explicit object dtype, so anndata's nullable-string guard no longer trips. This
# keeps any *other* nullable-string write in an aligner env from aborting a sweep.
# anndata reads it at import time; `conda run` passes the environment through, and
# files written this way read back correctly in aim_env.
export ANNDATA_ALLOW_WRITE_NULLABLE_STRINGS="${ANNDATA_ALLOW_WRITE_NULLABLE_STRINGS:-1}"

[[ -f "$PAIRS_CSV" ]] || die "pairs csv not found: $PAIRS_CSV"
[[ -d "$SC_DIR" && -d "$ST_DIR" ]] || die "expected $SC_DIR and $ST_DIR"
mkdir -p "$LOG_ROOT" "$WORK_DIR"

RUN_ID="$(date '+%Y%m%d-%H%M%S')-$$"
COUNTER_DIR="$LOG_ROOT/.counters/$RUN_ID"     # per-worker tallies for the summary

# A second instance for the same output + shard would repeat this one's sweeps and
# write into the same run roots (the split is per invocation - there is deliberately
# no coordination at runtime). One startup check against a recorded PID; a leftover
# file from a killed run is detected as stale and replaced.
RUNFILE="$LOG_ROOT/.running.shard${SHARD_INDEX}of${SHARD_COUNT}.pid"
if [[ "$DRY_RUN" != "1" ]]; then
  other=""
  [[ -f "$RUNFILE" ]] && other="$(tr -dc '0-9' < "$RUNFILE")"
  if [[ -n "$other" ]] && kill -0 "$other" 2>/dev/null; then
    if [[ "$FORCE" != "1" ]]; then
      die "another run for $OUT_ROOT (shard $SHARD_INDEX/$SHARD_COUNT) is active as pid $other.
       Two runs duplicate each other's sweeps into the same run roots. Stop it
       first ('kill $other'), or pass FORCE=1 if you are sure they are disjoint
       (different SHARD_INDEX, different OUT_ROOT)."
    fi
    log "FORCE=1: starting although pid $other is still running - make sure the two are disjoint"
  else
    [[ -n "$other" ]] && log "replacing stale $RUNFILE (pid $other is gone)"
    echo "$$" > "$RUNFILE"
    trap 'rm -f "$RUNFILE"' EXIT
  fi
fi

# --- group enumeration + static split ---------------------------------------
# Emits one line per (pair, start mode), prefixed with the worker it is assigned
# to:  worker \t PairID \t scName \t stName \t start_tag \t ann_column
# ("-" = Leiden over-clustering). Annotation columns come from scRNA/index.csv
# (CellTypeKey0/1/2); the built-in table is the fallback if that file is absent.
# Duplicate groups are dropped - a repeated pairs.csv row would otherwise become
# two work items pointing at one run root. Done in python so quoted CSV fields and
# CRLF line ends are handled correctly.
enumerate_and_split() {
  LEIDEN_RES="$LEIDEN_RES" START_MODES="$START_MODES" NO_ANN="$NO_ANN" ORDER="$ORDER" \
  WORKERS="$WORKERS" SHARD_COUNT="$SHARD_COUNT" SHARD_INDEX="$SHARD_INDEX" \
  python3 - "$PAIRS_CSV" "$SC_INDEX" "$ST_INDEX" <<'PY'
import csv, os, sys

pairs_csv, sc_index, st_index = sys.argv[1], sys.argv[2], sys.argv[3]
modes = os.environ["START_MODES"].split()
leiden_tag = "leiden_res" + os.environ["LEIDEN_RES"]
no_ann = os.environ["NO_ANN"]
order = os.environ.get("ORDER", "lpt")
workers = int(os.environ["WORKERS"])
shard_count, shard_index = int(os.environ["SHARD_COUNT"]), int(os.environ["SHARD_INDEX"])

# Fallback = the CellTypeKey0/1/2 columns of 01_Datasets/scRNA/index.csv.
FALLBACK = {
    "01_mouse-ssp":           ["cellType", "cellTypeMinor"],
    "02_colorectal-cancer":   ["cellType"],
    "03_mouse-liver":         ["celltype", "cell_ontology_class"],
    "04_mouse-kidney":        ["celltype", "free_annotation"],
    "05_mouse-embryo":        ["celltype"],
    "06_mouse-hypothalamic":  ["celltype"],
    "07_mouse-retina":        ["Cell.type", "cell_ontology_class"],
    "08_mouse-mop":           ["celltype"],
    "09_mouse-mpfc":          ["celltype"],
    "10_mouse-visual-cortex": ["celltype", "subclass"],
    "11_human-breast-cancer": ["cellType"],
    "12_human-mtg":           ["cellType"],
}


def index_column(path, name_col, value_col, cast):
    """{Name: value} from an index.csv column; empty when file/column is absent."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            name = (row.get(name_col) or "").strip()
            try:
                out[name] = cast(row.get(value_col))
            except (TypeError, ValueError):
                continue
    return out


ann_keys = {}
if os.path.exists(sc_index):
    with open(sc_index, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            name = (row.get("Name") or "").strip()
            keys = [(row.get(f"CellTypeKey{i}") or "").strip() for i in range(3)]
            keys = [k for k in keys if k]
            if name and keys:
                ann_keys[name] = keys
else:
    print(f"WARNING: {sc_index} missing - using the built-in annotation table",
          file=sys.stderr)

n_cells = index_column(sc_index, "Name", "CellCount", lambda v: int(float(str(v).strip())))
n_spots = index_column(st_index, "Name", "SpotCount", lambda v: int(float(str(v).strip())))

groups, seen, warned = [], set(), set()
with open(pairs_csv, newline="", encoding="utf-8-sig") as fh:
    for row in csv.DictReader(fh):
        pid = int(str(row["PairID"]).strip())
        sc, st = row["scName"].strip(), row["stName"].strip()
        if shard_count > 1 and pid % shard_count != shard_index:
            continue
        # Cost proxy for the split: a sweep scales with the reference and the
        # slice. Annotation modes get a discount - far fewer start clusters, hence
        # far fewer K levels, than Leiden at resolution 3.
        weight = n_cells.get(sc, 0) * max(n_spots.get(st, 0), 1)
        wanted = []
        if "leiden" in modes:
            wanted.append((leiden_tag, no_ann, weight))
        if "annotation" in modes:
            keys = ann_keys.get(sc, FALLBACK.get(sc, []))
            if not keys and sc not in warned:
                warned.add(sc)
                print(f"WARNING: no cell-type annotation known for {sc}", file=sys.stderr)
            for key in dict.fromkeys(keys):          # de-duplicate repeated keys
                wanted.append((f"ann_{key}", key, weight // 4))
        for tag, ann, w in wanted:
            key = (pid, sc, st, tag)
            if key in seen:                          # repeated pairs.csv row
                print(f"WARNING: duplicate group {key} ignored", file=sys.stderr)
                continue
            seen.add(key)
            groups.append((w, pid, sc, st, tag, ann))

if not groups:
    print("ERROR: no groups enumerated", file=sys.stderr)
    raise SystemExit(1)

# Static split, decided here and never revisited: greedy longest-processing-time
# first (each group to the currently least-loaded worker) if we have cost
# estimates, plain round-robin otherwise.
loads = [0] * workers
assigned = []
if order == "lpt":
    groups.sort(key=lambda g: -g[0])
    for g in groups:
        w = min(range(workers), key=lambda i: (loads[i], i))
        loads[w] += g[0]
        assigned.append((w + 1, g))
else:
    for i, g in enumerate(groups):
        w = i % workers
        loads[w] += g[0]
        assigned.append((w + 1, g))

for w, (_, pid, sc, st, tag, ann) in assigned:
    print(f"{w}\t{pid}\t{sc}\t{st}\t{tag}\t{ann}")

total = sum(loads) or 1
for i in range(workers):
    n = sum(1 for w, _ in assigned if w == i + 1)
    print(f"  worker {i + 1:>2}: {n:>4} groups, {100.0 * loads[i] / total:5.1f}% of the "
          f"estimated load", file=sys.stderr)
PY
}

# --- preflight: do the annotation columns exist in the h5ads? ---------------
# Reads only the /obs group's member names via h5py - never the count matrix.
preflight_annotations() {
  SC_DIR="$SC_DIR" NO_ANN="$NO_ANN" python3 - "$GROUPS_FILE" <<'PY'
import os, sys, h5py

groups_file = sys.argv[1]
sc_dir, no_ann = os.environ["SC_DIR"], os.environ["NO_ANN"]

wanted: dict[str, set[str]] = {}
with open(groups_file, encoding="utf-8") as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 6 or parts[5] == no_ann:
            continue
        wanted.setdefault(parts[2], set()).add(parts[5])

problems = 0
for sc, keys in sorted(wanted.items()):
    path = os.path.join(sc_dir, f"{sc}.h5ad")
    if not os.path.exists(path):
        print(f"  MISSING FILE  {path}")
        problems += 1
        continue
    try:
        with h5py.File(path, "r") as fh:
            cols = set(fh["obs"].keys())
    except Exception as exc:  # unreadable / not an h5ad
        print(f"  UNREADABLE    {path}: {exc}")
        problems += 1
        continue
    for key in sorted(keys):
        ok = key in cols
        problems += 0 if ok else 1
        print(f"  {'OK  ' if ok else 'FAIL'}  {sc}: obs[{key!r}]")
        if not ok:
            print(f"        available obs columns: {sorted(cols)}")
print(f"  -> {problems} annotation problem(s)")
PY
}

pair_dir_name() { printf '%03d_%s__%s' "$1" "$2" "$3"; }

# The reference scaffold (start clusters + aggregates + UMAPs) does not depend on
# the tree method, but it is cached per run root - and a run root has to carry the
# tree method. Hardlink it across sibling tree cells instead of recomputing it.
# Safe: a group belongs to exactly one worker, in every stage.
link_scaffold() {
  local cell="$1" pair_root="$2" prefix="$3" src
  [[ "$LINK_SCAFFOLD" == "1" ]] || return 0
  [[ -f "$cell/reference_scaffold.h5ad" ]] && return 0
  for src in "$pair_root"/*/reference_scaffold.h5ad; do
    [[ -f "$src" && -f "${src%.h5ad}.meta.json" ]] || continue
    mkdir -p "$cell"
    if ln "$src" "$cell/reference_scaffold.h5ad" 2>/dev/null &&
       ln "${src%.h5ad}.meta.json" "$cell/reference_scaffold.meta.json" 2>/dev/null; then
      log "$prefix scaffold hardlinked from $(dirname "$src")"
    else
      rm -f "$cell/reference_scaffold.h5ad" "$cell/reference_scaffold.meta.json"
    fi
    return 0
  done
}

ledger_row() {   # status seconds mapping start_tag tree pair_id sc st run_root
  printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
    "$(date '+%F %T')" "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" >> "$LEDGER"
}

# --- one sweep --------------------------------------------------------------
# Sets SWEEP_STATUS to ok | skipped | failed for the caller's counters.
run_sweep() {
  local mapping="$1" pid="$2" sc="$3" st="$4" tag="$5" ann="$6" tree="$7" prefix="$8"
  [[ "$ann" == "$NO_ANN" ]] && ann=""

  local pair_root cell run_root
  pair_root="$OUT_ROOT/$tag/$(pair_dir_name "$pid" "$sc" "$st")"
  cell="$pair_root/$tree"
  run_root="$cell/$mapping"

  if [[ -f "$run_root/k_comparison.csv" ]]; then
    SWEEP_STATUS=skipped
    ledger_row skipped 0 "$mapping" "$tag" "$tree" "$pid" "$sc" "$st" "$run_root"
    return 0
  fi

  local sc_path="$SC_DIR/$sc.h5ad" st_path="$ST_DIR/$st.h5ad"
  if [[ ! -f "$sc_path" || ! -f "$st_path" ]]; then
    log "$prefix MISSING INPUT - skipping ($sc_path / $st_path)"
    SWEEP_STATUS=failed
    ledger_row missing_input 0 "$mapping" "$tag" "$tree" "$pid" "$sc" "$st" "$run_root"
    return 0
  fi

  log "$prefix $mapping | $tag | $tree | pair $pid ($sc x $st)"
  mkdir -p "$cell"
  link_scaffold "$cell" "$pair_root" "$prefix"

  local args=(run
    --scdata "$sc_path"
    --stdata "$st_path"
    --output_dir "$cell"
    --mapping "$mapping"
    --agglo_tree_method "$tree"
    --leiden_resolution "$LEIDEN_RES"
  )
  [[ -n "$ann"    ]] && args+=(--start_from_annotation "$ann")
  [[ -n "$K_MIN"  ]] && args+=(--k_min "$K_MIN")
  [[ -n "$K_MAX"  ]] && args+=(--k_max "$K_MAX")
  [[ -n "$K_STEP" ]] && args+=(--k_step "$K_STEP")

  local job_log_dir job_log t0 dt rc
  job_log_dir="$LOG_ROOT/$tag/$(pair_dir_name "$pid" "$sc" "$st")/$tree"
  mkdir -p "$job_log_dir"
  job_log="$job_log_dir/$mapping.log"

  t0="$(date +%s)"
  if [[ -n "$JOB_TIMEOUT" ]]; then
    timeout "$JOB_TIMEOUT" aim "${args[@]}" > "$job_log" 2>&1
  else
    aim "${args[@]}" > "$job_log" 2>&1
  fi
  rc=$?
  dt=$(( $(date +%s) - t0 ))

  if (( rc == 0 )); then
    SWEEP_STATUS=ok
    log "$prefix ok in ${dt}s -> $run_root"
    ledger_row ok "$dt" "$mapping" "$tag" "$tree" "$pid" "$sc" "$st" "$run_root"
  else
    SWEEP_STATUS=failed
    log "$prefix FAILED (rc=$rc) after ${dt}s - see $job_log"
    tail -n 5 "$job_log" | sed "s|^|$prefix   \| |"
    ledger_row "failed_rc$rc" "$dt" "$mapping" "$tag" "$tree" "$pid" "$sc" "$st" "$run_root"
  fi
}

# --- one worker -------------------------------------------------------------
# Walks its own work file top to bottom. No shared state with the other workers.
run_worker() {
  local mapping="$1" wid="$2" work_file="$3" n_own="$4"
  local idx=0 done_n=0 skipped_n=0 failed_n=0
  local pid sc st tag ann tree prefix

  while IFS=$'\t' read -r pid sc st tag ann; do
    [[ -n "${pid:-}" ]] || continue
    idx=$(( idx + 1 ))
    for tree in $TREES; do
      prefix="[w$wid $mapping $idx/$n_own]"
      SWEEP_STATUS=failed
      run_sweep "$mapping" "$pid" "$sc" "$st" "$tag" "$ann" "$tree" "$prefix"
      case "$SWEEP_STATUS" in
        ok)      done_n=$(( done_n + 1 )) ;;
        skipped) skipped_n=$(( skipped_n + 1 )) ;;
        *)       failed_n=$(( failed_n + 1 )) ;;
      esac
    done
    # A scaffold cache is a full prepared reference; after the last mapper stage
    # of this group nothing needs it any more.
    if [[ "$PRUNE_SCAFFOLD" == "1" && "$mapping" == "$LAST_MAPPING" ]]; then
      for tree in $TREES; do
        rm -f "$OUT_ROOT/$tag/$(pair_dir_name "$pid" "$sc" "$st")/$tree/reference_scaffold.h5ad" \
              "$OUT_ROOT/$tag/$(pair_dir_name "$pid" "$sc" "$st")/$tree/reference_scaffold.meta.json"
      done
    fi
  done < "$work_file"

  echo "$done_n $skipped_n $failed_n" > "$COUNTER_DIR/$mapping.$wid"
}

# --- build the split --------------------------------------------------------
log "Enumerating groups and splitting them across $WORKERS worker(s)..."
enumerate_and_split > "$GROUPS_FILE.tmp" || die "group enumeration failed"
mv "$GROUPS_FILE.tmp" "$GROUPS_FILE"
[[ -s "$GROUPS_FILE" ]] || die "no groups enumerated from $PAIRS_CSV"

rm -f "$WORK_DIR"/worker_*.tsv
for (( wid = 1; wid <= WORKERS; wid++ )); do
  awk -F'\t' -v w="$wid" '$1 == w { print $2"\t"$3"\t"$4"\t"$5"\t"$6 }' \
    "$GROUPS_FILE" > "$WORK_DIR/worker_$wid.tsv"
done

N_GROUPS="$(wc -l < "$GROUPS_FILE" | tr -d ' ')"
N_ASSIGNED="$(cat "$WORK_DIR"/worker_*.tsv | wc -l | tr -d ' ')"
(( N_GROUPS == N_ASSIGNED )) || die "split lost groups: $N_GROUPS enumerated, $N_ASSIGNED assigned"

N_TREES="$(printf '%s\n' $TREES | wc -l | tr -d ' ')"
N_MAPPERS="$(printf '%s\n' $MAPPINGS | wc -l | tr -d ' ')"
TOTAL=$(( N_GROUPS * N_TREES * N_MAPPERS ))
LAST_MAPPING="$(printf '%s\n' $MAPPINGS | tail -1)"

K_DESC="full sweep (K = n_start_clusters .. 1, step 1)"
if [[ -n "$K_MIN$K_MAX$K_STEP" ]]; then
  K_DESC="${K_MIN:+--k_min $K_MIN }${K_MAX:+--k_max $K_MAX }${K_STEP:+--k_step $K_STEP}"
fi

cat <<EOF

================================ Experiment 1 ================================
  pairs csv    : $PAIRS_CSV
  data         : $SC_DIR | $ST_DIR
  output       : $OUT_ROOT
  mappers      : $MAPPINGS   (one stage each, in this order)
  trees        : $TREES
  start modes  : $START_MODES (leiden resolution $LEIDEN_RES)
  K range      : $K_DESC
  shard        : $SHARD_INDEX / $SHARD_COUNT   (by PairID)
  groups       : $N_GROUPS   (pairs x start modes, split $ORDER across the workers)
  sweeps       : $TOTAL   (groups x trees x mappers)
  workers      : $WORKERS   ($OMP_NUM_THREADS thread(s) each, $NPROC cores seen)
  work files   : $WORK_DIR/worker_<n>.tsv
  run id       : $RUN_ID
==============================================================================

EOF

if (( WORKERS > 1 )); then
  log "NOTE: each worker loads its own sc + ST AnnData - size WORKERS by RAM, not by cores."
fi

if [[ "$DRY_RUN" == "1" ]]; then
  log "Annotation preflight (reads only obs column names, never the matrix):"
  preflight_annotations
  log "Groups per worker:"
  for (( wid = 1; wid <= WORKERS; wid++ )); do
    printf '  worker %2d: %4d groups, first: %s\n' "$wid" \
      "$(wc -l < "$WORK_DIR/worker_$wid.tsv" | tr -d ' ')" \
      "$(head -1 "$WORK_DIR/worker_$wid.tsv" | awk -F'\t' '{print $4" "$1" "$2" x "$3}')"
  done
  log "DRY_RUN=1 - nothing was executed."
  exit 0
fi

[[ -f "$LEDGER" ]] || echo "timestamp,status,seconds,mapping,start_tag,tree,pair_id,sc,st,run_root" > "$LEDGER"
mkdir -p "$COUNTER_DIR"

# Job control, so every worker becomes its own process group and Ctrl-C / SIGTERM
# can take down the worker *and* the aim process it is currently running.
set -m
declare -a WORKER_PIDS=()

stop_workers() {
  trap - INT TERM
  log "interrupted - terminating $WORKERS worker(s)"
  local p
  for p in "${WORKER_PIDS[@]:-}"; do
    [[ -n "$p" ]] || continue
    kill -TERM "-$p" 2>/dev/null || kill -TERM "$p" 2>/dev/null
  done
  wait
  log "stopped. Re-run the same command to resume."
  exit 130
}
trap stop_workers INT TERM

# --- run: one stage per mapper, each worker on its own groups ---------------
started_at="$(date +%s)"
grand_done=0 grand_skipped=0 grand_failed=0

for mapping in $MAPPINGS; do
  stage_t0="$(date +%s)"
  log "=== stage '$mapping': $N_GROUPS groups over $WORKERS worker(s)"

  WORKER_PIDS=()
  for (( wid = 1; wid <= WORKERS; wid++ )); do
    n_own="$(wc -l < "$WORK_DIR/worker_$wid.tsv" | tr -d ' ')"
    (( n_own > 0 )) || continue
    run_worker "$mapping" "$wid" "$WORK_DIR/worker_$wid.tsv" "$n_own" &
    WORKER_PIDS+=("$!")
  done
  wait

  s_done=0 s_skipped=0 s_failed=0
  for f in "$COUNTER_DIR/$mapping".*; do
    [[ -f "$f" ]] || continue
    read -r d s fl < "$f"
    s_done=$(( s_done + d )); s_skipped=$(( s_skipped + s )); s_failed=$(( s_failed + fl ))
  done
  grand_done=$(( grand_done + s_done ))
  grand_skipped=$(( grand_skipped + s_skipped ))
  grand_failed=$(( grand_failed + s_failed ))

  stage_dt=$(( $(date +%s) - stage_t0 ))
  log "=== stage '$mapping' done in $(( stage_dt / 3600 ))h $(( (stage_dt % 3600) / 60 ))m: $s_done ok, $s_skipped skipped, $s_failed failed"
done

elapsed=$(( $(date +%s) - started_at ))
cat <<EOF

============================== Experiment 1 done =============================
  sweeps      : $TOTAL
  completed   : $grand_done
  skipped     : $grand_skipped  (already had k_comparison.csv)
  failed      : $grand_failed
  workers     : $WORKERS
  wall clock  : $(( elapsed / 3600 ))h $(( (elapsed % 3600) / 60 ))m
  ledger      : $LEDGER
  logs        : $LOG_ROOT
==============================================================================
EOF

if (( grand_failed > 0 )); then
  log "Failing combinations (from the ledger):"
  awk -F, '$2 ~ /^(failed|missing_input)/ {printf "  %-14s %-16s %-22s %-8s pair %s\n",$2,$4,$5,$6,$7}' \
    "$LEDGER" | sort -u
  exit 1
fi