# Experiment 1 — on NHR@FAU TinyFat

Every pair of `Data/01_Datasets_Paper/pairs.csv` × every start mode × both linkages × all five
mappers, through the `spatialaim` CLI, with **wall time and peak memory measured per run** by the
runner (SpatialAIM itself is not touched). Sections 03 and 04 read the pulled summary.

## The grid

| | |
|---|---|
| start modes | `leiden_res3.0` + `ann_<key>` for every `CellTypeKey*` of the reference → 88 + 164 = **252 (pair, start) units** |
| stage A — TinyFat | one reference-scaffold build per unit (252) |
| stage B — TinyFat | {ward, average} × {nearest_centroid, wann, tacco, dot, tangram} per unit (2520); Tangram on CPU |
| K | full sweep, step 1; Leiden resolution 3.0 in both modes |

**Why two stages.** The reference scaffold (start clusters, shared-gene Leiden, aggregates, both
UMAPs) depends on the pair, the resolution and the start mode, not on the linkage or the mapper. Stage
A builds it once per unit (`nearest_centroid` at K = 1, everything but the scaffold discarded);
each stage-B run copies it into its own fresh output dir, so its time and memory are the sweep
alone. `run_one.py` checks every stage-B log for `Reusing cached reference scaffold`. The summary
carries both parts, so the full cost of a configuration is `wall_s + scaffold_build_s`.

**Measurement** (`run_one.py`, per run → `resources.json`):
- wall time, user and system CPU time;
- peak RSS summed over the whole process tree, sampled every 0.5 s, including the `conda run`
  aligner children; also the largest single process and the children alone;
- host, CPU model, GPU name, thread count, memory limit, job ids and the spatialaim version.

Every job gets 8 physical cores and 64 GB, and OMP/MKL/OpenBLAS/numba are pinned to 8 threads.
All mappers run on the same hardware (TinyFat `work`, AMD EPYC 7502), so their timings compare
directly. Tangram runs on CPU: GPU jobs on TinyGPU were rarely scheduled, and `run_tangram.py`
falls back to `cpu` when no GPU is visible. State that in the paper (Tangram is usually run on a
GPU).

**Caveat on the timings.** TinyFat `work` nodes are shared: cores are exclusive, memory bandwidth
is not.

**Jobs.** The grid is fixed up front in `units.tsv` (checked in, written by `grid.py`): line i
is unit i, and i is the Slurm array index of every job on that unit.

| array | one job = | jobs | runs per job |
|---|---|---|---|
| `01_scaffolds.tinyfat.sbatch` | one unit's scaffold | 252 | 1 |
| `02_map.tinyfat.sbatch <mapper>`, once per mapper | (unit, mapper): ward, then average | 5 × 252 | 2 |

The slowest run of the first round was a 7.3 h DOT sweep, so two fit the 24 h limit. Units are
ordered by PairID; units 0–2 are pair 1 (`01_mouse-ssp`), which the one-pair test runs.

**Resume.** Every run is idempotent: a `DONE` marker means skip, a `FAILED` marker means retry.
Anything killed by the time limit simply has no marker. Nothing is ever deleted except a run's own
scratch.

## Files

| File | Role |
|---|---|
| `00_setup_server.sh` | clone SpatialAIM-Repository at the pinned commit; create `spatialaim_env`, `tangram_env`, `tacco_env`, `dot_env` (+ DOTr); `check` mode runs smoke checks |
| `01_scaffolds.tinyfat.sbatch`, `02_map.tinyfat.sbatch` | the arrays; first argument `run` or `test` (the run root), then (stage B) the mapper |
| `03_status.py` | done / failed / missing per stage and mapper, cache misses, OOM suspects, and the `sbatch --array=…` lines that resubmit what is unfinished |
| `04_summarize.py` | → `summary/` (below) |
| `05a_run_k.py`, `05a_run_k.tinyfat.sbatch`, `05b_merge_k.py` | only for DOT on unit 248 (pair 87, breast cancer, Leiden start), which does not fit 24 h: one job per K, then merged into the normal run folder (step 5b) |
| *helpers, not run by hand* | |
| `config.sh` | shared settings, sourced by everything: `$WORK/spatialaim_exp1` layout, conda in `$WORK`, proxy, caches off `$HOME`, thread pinning |
| `units.tsv`, `grid.py` | the 252 units and the walltimes; `python grid.py array leiden\|annotation [--test]` prints an `--array` list; `python grid.py units` rewrites `units.tsv` (only if the datasets change) |
| `run_one.py` | the measured runner: `run_one.py A <unit>` or `run_one.py B <unit> <mapper>` |
| `run_job.sh` | body of every array task |

## On the cluster

```
$WORK/spatialaim_exp1/
├── code/SpatialAIM-Repository/     pinned clone (installed non-editable into spatialaim_env)
├── code/experiment1/               this folder; submit from here
│   └── logs/                       Slurm output, <job name>_<array job>_<unit>.out
├── data/01_Datasets_Paper/         read in place by both stages - never move or re-copy between A and B
├── envs/  conda_pkgs/  cache/
├── test/                           run root of the one-pair test
└── run/                            run root of the full grid
    ├── scaffolds/<PID>_<sc>__<st>/<start_tag>/    reference_scaffold.{h5ad,meta.json}, resources.json, run.log, DONE
    ├── results/<PID>_<sc>__<st>/<start_tag>/<tree>/<mapper>/
    │       config.yaml, k_comparison.csv, k_selection.json,
    │       k_levels.zip  (every k_<kkk>/ file except the per-spot ones: tree cut, analysis/data JSONs,
    │                      per-gene cossim; read with Python's zipfile, no need to unpack),
    │       run.log, resources.json, DONE | FAILED
    ├── logs/       ledger.tsv (one line per run)
    └── summary/
```

The scaffolds take roughly 1–2 GB each for the large references, so tens of GB in total; the
results are small. `$WORK` gives 1 TB.

## Step by step

Placeholders: `<user>` is your HPC account, and `<WORK>` is what `echo $WORK` prints on tinyx
(e.g. `/home/woody/<group>/<user>`). Local commands are PowerShell; cluster commands are bash on
tinyx.

### 0. SSH (once, local)

The SSH key is uploaded in the HPC portal. Add to `~/.ssh/config`:

```
Host csnhr
    HostName csnhr.nhr.fau.de
    User <user>
    IdentityFile ~/.ssh/id_ed25519_nhr_fau
    IdentitiesOnly yes
Host tinyx
    HostName tinyx.nhr.fau.de
    User <user>
    ProxyJump csnhr
    IdentityFile ~/.ssh/id_ed25519_nhr_fau
    IdentitiesOnly yes
```

### 1. Copy data and code (local)

`csnhr` is the transfer host.

```powershell
ssh tinyx "mkdir -p <WORK>/spatialaim_exp1/data <WORK>/spatialaim_exp1/code; echo `$WORK"
scp -r C:\Users\zi69hebi\Dev\10_Alignment\Data\01_Datasets_Paper csnhr:<WORK>/spatialaim_exp1/data/
scp -r C:\Users\zi69hebi\Dev\10_Alignment\SpatialAIM-Experiments\experiment1 csnhr:<WORK>/spatialaim_exp1/code/
```

After a later code change, copy `experiment1/` again; the data never.

### 2. Environments (tinyx, once, takes a while)

```bash
cd $WORK/spatialaim_exp1/code/experiment1
bash 00_setup_server.sh                     # clone @ pinned commit, all four envs, checks
bash 00_setup_server.sh check               # checks only, any time later
```

The check has to show `OK` for all three wrappers and for DOTr, and `spatialaim validate` must
show no errors.

### 3. One-pair end-to-end test

This runs pair 1 (`01_mouse-ssp`) with one start of each kind into its own root `test/`: unit 0
(Leiden start) and unit 1 (`ann_cellType`), so 2 scaffolds and 5 mappers × 2 trees × 2 units = 20
runs. That covers both code paths and both roles of `best_per_run.csv` (`TEST_UNITS` in `grid.py`). Every
command runs in `$WORK/spatialaim_exp1/code/experiment1`.

```bash
cd $WORK/spatialaim_exp1/code/experiment1
source config.sh && conda activate spatialaim_env
export EXP1_ROOT=$EXP1_BASE/test

sbatch.tinyfat --array=0-1 01_scaffolds.tinyfat.sbatch test
squeue.tinyfat -u $USER                     # wait until they are gone (a few minutes)
python 03_status.py                         # A: 2 done

for m in nearest_centroid wann tacco tangram; do
  sbatch.tinyfat -J exp1_B_$m --time=02:00:00 --array=0-1 02_map.tinyfat.sbatch test $m
done
sbatch.tinyfat -J exp1_B_dot --time=08:00:00 --array=0 02_map.tinyfat.sbatch test dot   # Leiden
sbatch.tinyfat -J exp1_B_dot --time=04:00:00 --array=1 02_map.tinyfat.sbatch test dot   # annotation
python 03_status.py --failed                # wait until 22 / 22 done (~1 h, mostly dot)

python 04_summarize.py
```

Check before going on:
- [ ] `03_status.py`: 22 / 22 done, 0 failed, **0 cache misses**.
- [ ] `summary/meta.json`: `missing: 0`.
- [ ] `summary/best_per_run.csv`: 20 rows, 10 `spatialaim` and 10 `baseline`.
- [ ] `summary/resources_stage_b.csv`: tacco, dot and tangram rows
      `peak_rss_children_mb` > 0 (the aligner process was sampled); tangram `wall_s` minutes, not
      seconds.
- [ ] `summary/resources_stage_a.csv`: `scaffold_s` well under the 30 min limit.
- [ ] `logs/*.out`: no tracebacks.

Pull the test (local) and draw it into a scratch vault, which checks the figure side:

```powershell
$W = "C:\Users\zi69hebi\Dev\10_Alignment\Data\02_Results\experiment1_test"; mkdir $W
scp -r "csnhr:<WORK>/spatialaim_exp1/test/summary" $W
$env:SPATIALAIM_EXP1 = $W; $env:SPATIALAIM_VAULT = "$env:TEMP\exp1_vault"
python sections/03_reconstruction/panels.py; python sections/04_biology/panels.py --mapper nearest_centroid
```

These figures are sparse because they have one reference only; the point is that the pipeline
runs.

### 4. Full grid

Walltimes (`JOB_TIME`, `DOT_TIME` in `grid.py`): 2 h for nearest_centroid, wann, tacco and tangram (first
round's slowest job: 0.7, 0.7, 1.3 and 0.9 h), so they can also use TinyFat's short-job nodes. DOT is split by start kind
(annotation starts sweep few K and are fast) and gets realistic limits instead of 24 h for all
252 jobs: 4 h for annotation starts (first round: p99 2.6 h, max 5.1 h) and 8 h for Leiden starts
(p99 6.6 h, max 8.2 h). The few jobs that run over get resubmitted with `--time=24:00:00` (step 5);
a tree that finished before has `DONE` and is skipped.

`%20` caps how many jobs of one array run at once (6 stage-B arrays → at most 120 jobs); raise or
lower it as the queue allows.

```bash
export EXP1_ROOT=$EXP1_BASE/run
sbatch.tinyfat --array=0-251%20 01_scaffolds.tinyfat.sbatch run
python 03_status.py                         # wait until A: 252 done

for m in nearest_centroid wann tacco tangram; do
  sbatch.tinyfat -J exp1_B_$m --time=02:00:00 --array=0-251%20 02_map.tinyfat.sbatch run $m
done
sbatch.tinyfat -J exp1_B_dot --time=04:00:00 --array=$(python grid.py array annotation)%20 02_map.tinyfat.sbatch run dot
sbatch.tinyfat -J exp1_B_dot --time=08:00:00 --array=$(python grid.py array leiden)%20 02_map.tinyfat.sbatch run dot
```

Submit stage B only once stage A is complete: a stage-B job whose scaffold is missing reports
`no_scaffold` and exits. To try one job per array first, submit `--array=0`, check its
runs, then the full range: done runs are skipped.

Watch the run with `squeue.tinyfat -u $USER`, `python 03_status.py`
and `tail -f $EXP1_ROOT/logs/ledger.tsv`.

### 5. Leftovers

```bash
python 03_status.py --failed                # OOM? / exit <n> with the log path, then one
                                            # ready sbatch line per array with the unfinished units
```

Copy those lines; a resubmitted job skips its runs that already have `DONE`. Extra options go
before the script name and override its `#SBATCH` lines, e.g. `--mem=128G` after an OOM or
`--time=…`.

### 5b. The one run that does not fit 24 h

DOT on unit 248 (pair 87, `09_human-breast-cancer` against the whole 167,780-cell Xenium
section, Leiden start) needs more than 24 h for one tree. It runs as **one array task per K**
instead: each task is `spatialaim run --k_min K --k_max K` on the unit's scaffold, measured like
any run, and `05b_merge_k.py` then scores the K levels together with spatialaim's own cross-K
comparison and writes the ordinary run folder (`k_comparison.csv`, `k_selection.json`,
`k_levels.zip`, `resources.json`, `DONE`), so `03_status.py` and `04_summarize.py` see a normal
run. Per tree, only for a tree without `DONE`:

```bash
# which trees already finished (a listed DONE = that tree needs nothing more):
ls $EXP1_BASE/run/results/087_09_human-breast-cancer__09_human-breast-cancer/leiden_res3.0/*/dot/DONE
sbatch.tinyfat --array=$(python 05a_run_k.py ward --levels)%20 05a_run_k.tinyfat.sbatch run ward
python 05b_merge_k.py ward                  # once all K are done; else prints the --array to resubmit
```

(and the same with `average`). Each K keeps its files and log under
`$EXP1_BASE/run/split_k/<run_id>/k_<kkk>/`. In the summary this run's `wall_s` is the **sum over its K
jobs** and its peak memory the largest of them (`resources.json` carries `split_by_k: true` and a
`per_k` list with each job's time, memory, mapping and analysis seconds). Every K job loads the scaffold and the slice on its own, so the sum is a little
above what one sweep would take -- say so where the timing is reported.

### 6. Summary and pull

```bash
python 04_summarize.py                      # expects missing: 0
```

```powershell
$D = "C:\Users\zi69hebi\Dev\10_Alignment\Data\02_Results\experiment1"; mkdir $D
scp -r "csnhr:<WORK>/spatialaim_exp1/run/summary" $D
```

### 7. Figures (local)

```powershell
python sections/03_reconstruction/panels.py           # panels into the vault (safe)
foreach ($m in "nearest_centroid", "wann", "tacco", "dot", "tangram") { python sections/04_biology/panels.py --mapper $m }
python sections/03_reconstruction/make_excalidraw.py  # OVERWRITES the hand-edited drawing
python sections/04_biology/make_excalidraw.py --mapper <m>  # idem, per mapper
```

Step `panels.py` alone refreshes the images the existing drawings embed. Run `make_excalidraw.py`
only if you want the layout regenerated.

## `summary/`

| File | Content |
|---|---|
| `best_per_run.csv` | one row per stage-B run root at the K that matters: `spatialaim` (Leiden, K = `best.overall`) or `baseline` (annotation, K = T). Only what sections 03/04 read: `pair_id`, `reference`, `pair_dir`, `role`, `ann_key`, `tree`, `mapper`, `k`, `score_overall` and the six plotted metrics |
| `all_k_per_run.csv.gz` | every `k_comparison.csv` row with the run's keys, for anything beyond the figures |
| `resources_stage_b.csv` | per stage-B run: `wall_s`, `cpu_user_s`, `cpu_sys_s`, `peak_rss_tree_mb`, `peak_rss_proc_mb`, `peak_rss_children_mb`, host/CPU/GPU, plus `scaffold_build_s`, `scaffold_peak_rss_tree_mb` and `wall_with_scaffold_s` |
| `resources_stage_a.csv` | per stage-A build: the same measurements plus `scaffold_s` (log-timed, to the cache write) and `scaffold_mb` |
| `missing.csv` | grid runs without `DONE` |
| `meta.json` | counts (grid runs, done, missing), spatialaim version(s) |
