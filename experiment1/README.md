# Experiment 1 on the server — install + run everything

The prerequisite sweep for the whole paper: every pair × start mode × linkage ×
mapper, one `aim run` each, ~2630 run roots. Sections 03 and 04 read nothing but
the metrics it collects; section 02 has its own small runs instead.

These scripts are self-contained — they install `spatial-aim` from TestPyPI on
the server and import nothing from this repository, so they are copied over as
they are.

Two scripts, no repository needed on the server:

| File | What it does |
|------|--------------|
| `01_setup_envs.sh` | creates all four envs — `aim_env` (with `spatial-aim` from **TestPyPI**) plus `tangram_env`, `tacco_env`, `dot_env` from the env specs shipped inside the wheel — then verifies each aligner env can run its wrapper |
| `02_run_experiment1.sh` | enumerates the whole grid and runs it across `WORKERS` shell workers; resumable, logged |
| `04_write_k_selection.py` | recomputes the missing `k_selection.json` per run root from `k_comparison.csv` alone, and reports truncated tables |
| `05_collect_metrics.sh` | tars the three small files per run root (~50 MB) for transfer back |

> `03_progress.sh`, the read-only progress report referred to below, was never
> copied back from the server and is not in this repository.

Once the archives are extracted into `Data/02_Results/experiment1`, everything
downstream goes through `aim_experiments.exp1`, which reduces each run root to
two rows (AIM's proposal, the annotation baseline) and caches the table in
`outputs/experiment1/runs.csv`.

## 0. Copy them over

```bash
scp -r aim_experiment1 <server>:/data/bionets/zi69hebi/Alignment/
```

## 1. Install

```bash
bash /data/bionets/zi69hebi/Alignment/aim_experiment1/01_setup_envs.sh
```

* **Nothing is ever deleted.** The script does not remove conda envs and never touches a results
  directory, so the finished `nearest_centroid` / `wann` sweeps stay where they are. An env that
  already exists is reused as-is — `conda env remove -n <env>` yourself to force a rebuild, or use
  `ONLY=tangram_env` to build a single env.
* `aim_env`: dependencies come from real PyPI, then `spatial-aim` from TestPyPI with `--no-deps`, so a
  TestPyPI mirror of a dependency can never be picked up by accident. Pin with
  `AIM_VERSION=0.1.4 bash 01_setup_envs.sh`.
* The aligner envs get **only their own aligner** — no `spatial-aim`. AIM runs their wrapper scripts
  by path (`conda run -n tangram_env python …/run_tangram.py …`), and those wrappers import nothing
  from `aim`. `dot_env` additionally gets the DOT R package from GitHub (`SKIP_DOT_R=1` to skip).
* pandas is pinned `<3` in `aim_env`: with pandas' new string default, AIM's own h5ad writes
  (`src/aim/io.py`) would hit anndata's nullable-string guard the way the aligner wrappers did. The
  wrappers are fixed; `io.py` is not yet.
* The last step prints `available: [...]` — all of `tangram`, `tacco`, `dot` must be listed — and then
  checks that each aligner env can actually execute its wrapper script. If conda is not on `PATH`
  inside the run shell, `export CONDA_EXE=...`.

## 2. Preflight (no compute)

```bash
DRY_RUN=1 bash /data/bionets/zi69hebi/Alignment/aim_experiment1/02_run_experiment1.sh
```

Prints the grid size, the start modes per reference, and checks that every annotation column really
exists in the corresponding `scRNA/*.h5ad` (reads only the `obs` column names, never the matrix).

## 3. Run in tmux

```bash
tmux new -s aim1 -d 'WORKERS=8 bash /data/bionets/zi69hebi/Alignment/aim_experiment1/02_run_experiment1.sh'
tmux attach -t aim1
```

Console lines are prefixed `[w<worker> <stage> g<group>/<total>]`. Progress:

```bash
tail -f /data/bionets/zi69hebi/Alignment/02_Results/experiment1/_logs/ledger.csv
```

Stop it: `Ctrl-C` in the session (the script terminates its workers and the `aim` process each one is
running, then tells you to re-run to resume). If the parent gets killed hard instead, workers survive
as orphans — `pkill -f 'aim run --scdata'` clears them. Either way nothing is corrupted: re-running
skips every sweep that already has a `k_comparison.csv`.

### How the work is split

All parallelism is in the shell script; AIM itself runs exactly as it does single-threaded. The split
is decided **before anything starts** and never revisited — there is no claiming, locking or other
coordination while the run is going.

* **Unit of work ("group") = one (pair, start mode)**, carrying both trees. A group maps 1:1 onto the
  directory `<start_tag>/<PairID>_<sc>__<st>/`, so two different groups can never touch the same
  directory. That is what makes the reference-scaffold cache safe: built once per group, then
  hardlinked from `ward` into `average`.
* **The group list is deduplicated** (a repeated `pairs.csv` row would otherwise become two work items
  pointing at one run root — it is reported as `WARNING: duplicate group … ignored`) and then **dealt
  out to the workers up front**, greedy longest-processing-time-first on
  `CellCount × SpotCount` from the index CSVs (annotation modes discounted, since they have far fewer K
  levels). Each worker gets `_logs/work/worker_<n>.tsv` and reads only that file, in every stage.
  The startup banner prints each worker's share, e.g. with `WORKERS=8`: one worker gets the single
  huge `11_human-breast-cancer` full-section group (≈18 % of the estimated load) and the other seven
  get 44–47 groups each (≈11.7 % apiece).
* **One stage per mapper**, in `MAPPINGS` order, with a barrier between stages: the two cheap
  in-process mappers finish across the entire grid first, and every later stage reads the scaffold
  cache the first stage wrote rather than building its own.
* **Only one run at a time per output + shard.** The script records its PID in
  `_logs/.running.shard<i>of<n>.pid` and refuses to start while that process is alive — a second
  invocation would otherwise repeat the first one's sweeps into the same run roots. A leftover file
  from a killed run is detected as stale and replaced. `FORCE=1` overrides it (use it for two shards
  on the same host).

`WORKERS=N` also splits the thread budget: `OMP_NUM_THREADS` defaults to `nproc / WORKERS`. Size `N`
by **RAM**, not by cores — every worker holds its own scRNA + ST AnnData.

For several machines, combine with the shard filter (disjoint by `PairID`, so shards never touch the
same group):

```bash
SHARD_COUNT=3 SHARD_INDEX=0 WORKERS=8 bash 02_run_experiment1.sh   # on host 1, etc.
```

## The grid

* **Pairs**: every row of `pairs_only_raw_counts.csv` (120 pairs; 10 of the 12 references —
  `06_mouse-hypothalamic` and `08_mouse-mop` do not appear in that subset).
* **Start clusters**: Leiden over-clustering (`--leiden_resolution 3.0`) **plus** every cell-type
  annotation the pair's reference has, taken from `scRNA/index.csv` `CellTypeKey0/1/2`
  (2 keys for 01, 03, 04, 07, 10; 1 key for 02, 05, 09, 11, 12) → **321 groups**.
* **Tree**: `ward`, `average`.
* **Mapper**: `nearest_centroid`, `wann`, `tangram`, `tacco`, `dot`.

**321 × 2 × 5 = 3210 sweeps**, each a full K sweep (K = number of start clusters down to 1).

Batch mode (`aim run --pairs_csv`) is deliberately not used: `--start_from_annotation` differs per
scRNA reference, and every start mode / tree combination needs its own output root.

## Output layout

```
02_Results/experiment1/
├── leiden_res3.0/<PairID>_<scName>__<stName>/{ward,average}/<mapper>/k_<kkk>/…
├── ann_<obs_column>/<PairID>_<scName>__<stName>/{ward,average}/<mapper>/k_<kkk>/…
└── _logs/
    ├── ledger.csv                        # one row per sweep: status, seconds, run root
    ├── work/groups_assigned.tsv          # the whole split, worker id in column 1
    ├── work/worker_<n>.tsv               # what worker n runs, in order
    ├── .running.shard0of1.pid            # PID of the active run (removed on exit)
    └── <start_tag>/<pair>/<tree>/<mapper>.log
```

A run root is named after the mapper alone, so the start mode and the tree method have to be
directory levels above it — that is exactly what `<start_tag>/…/<tree>/` provides. Point `aim gui`
at any `…/<tree>/` folder and it sees all five mappers side by side.

## Knobs

| Env var | Default | Notes |
|---------|---------|-------|
| `WORKERS` | `1` | parallel sweeps; `1` = plain sequential |
| `MAPPINGS` | `nearest_centroid wann tangram tacco dot` | also the stage order |
| `TREES` | `ward average` | |
| `START_MODES` | `leiden annotation` | drop one to run only the other |
| `LEIDEN_RES` | `3.0` | part of the `leiden_res*` folder name |
| `K_MIN` / `K_MAX` / `K_STEP` | unset → full sweep, step 1 | see the runtime note below |
| `ORDER` | `lpt` | weighted greedy split; `file` = `pairs.csv` order, dealt round-robin |
| `FORCE` | `0` | start although another run for this output + shard is alive |
| `OMP_NUM_THREADS` | `nproc / WORKERS` | also sets OPENBLAS/MKL/NUMBA |
| `JOB_TIMEOUT` | unset | e.g. `12h`, kills a single hanging sweep and moves on |
| `SHARD_COUNT` / `SHARD_INDEX` | `1` / `0` | across machines, by `PairID` |
| `PRUNE_SCAFFOLD` | `0` | delete `reference_scaffold.h5ad` after a group's last mapper stage |
| `LINK_SCAFFOLD` | `1` | hardlink the ward scaffold cache into the average cell |
| `DRY_RUN` | `0` | enumerate + check annotations only |
| `CACHE_ROOT` / `CONFIG_ROOT` | `<DATA_ROOT>/../.cache` / `.config` | matplotlib/numba/XDG caches, kept off `$HOME` |
| `OUT_ROOT`, `DATA_ROOT`, `PAIRS_CSV`, `AIM_ENV` | server defaults | |

## Fixed since the first run

**Nothing writes to `$HOME` any more.** The 2026-08-12 tangram run lost 582 of 634 sweeps to this:
the home directory went unwritable at 22:54 and from that minute on every `conda run -n tangram_env`
died before Python started (`An unexpected error has occurred`, because conda writes into
`~/.conda`), surfacing as `RuntimeError: tangram aligner failed (exit 1)` at the sweep's first K.
matplotlib hit the same wall one line earlier and only warned — `Permission denied:
'/home/…/.config/matplotlib'` is the fingerprint in the logs. Both scripts now export
`MPLCONFIGDIR` / `XDG_CACHE_HOME` / `XDG_CONFIG_HOME` / `NUMBA_CACHE_DIR` under `CACHE_ROOT` and pin
`CONDA_PKGS_DIRS` / `CONDA_ENVS_DIRS` to `$CONDA_BASE`, and the runner refuses to start if any of
those is unwritable, rather than burning through the grid at ~50 s per failed sweep.

**Nullable-string writes.** Every `tangram`/`tacco`/`dot` sweep used to die at its first K with
`anndata.settings.allow_write_nullable_strings is False`, because the wrappers built their
`mapping_prob.h5ad` index from Python strings and pandas' new string default made it a nullable
`str` extension array. The wrappers now pass `pd.Index(…, dtype=object)` explicitly. The runner still
exports `ANNDATA_ALLOW_WRITE_NULLABLE_STRINGS=1` as a net for any other write in those envs.

**No `spatial-aim` in the aligner envs.** The wrappers were launched as `python -m
aim.reference_aligners.run_tangram`, which required the package inside `tangram_env`. They are now
standalone scripts executed by path, importing nothing from `aim`, so the aligner envs hold only
their own aligner. The persistent per-sweep worker (`--server` mode, `AlignerWorker`) is gone with
it: each K is one `conda run`.

## Two things worth deciding before you start it

**Runtime.** In Leiden mode at resolution 3.0 a sweep has as many K levels as there are start
clusters (often > 100), each level running its own analysis — and for `tangram`/`tacco`/`dot` its own
external alignment. 1926 of the 3210 sweeps are external-aligner sweeps, and 71 of the 120 pairs are
`07_mouse-retina` slices. Since the persistent worker is gone, every one of those K levels now also
pays a fresh `conda run` (env activation + the aligner's imports + a re-read of the prepared pair) —
roughly 10–20 s each on top of the alignment itself. Even with 8 workers, full-grid full-K is a long
run. Options: coarsen the external-aligner stages (`K_STEP=5`, or `K_MAX` to cap the top), raise
`WORKERS`, or add machines via `SHARD_COUNT`.

Resume is per sweep, so none of this is lost work: anything with a `k_comparison.csv` is skipped, and
the finished `nearest_centroid` / `wann` sweeps from the first run are simply skipped on the next
start. Restrict the run to what is missing with `MAPPINGS="tangram tacco dot"`.

**Disk.** The reference-scaffold cache is one prepared reference h5ad per (pair, start mode, tree)
cell. Hardlinking ward → average halves that, leaving 321 copies. If the reference h5ads are large,
set `PRUNE_SCAFFOLD=1` — the GUI rebuilds a scaffold on demand when you browse a cell later.