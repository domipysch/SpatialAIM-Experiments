# Experiment 1 — on NHR@FAU TinyFat

Every pair × start × tree × mapper through `spatialaim run`, with time and memory measured per run,
per phase and per K.

## Grid

| | |
|---|---|
| units | (pair, start): `leiden` + `annotation_<key>` per `CellTypeKey*` → **193** for 94 pairs (`units.tsv`) |
| scaffold jobs | one reference scaffold per pair (94), reused by all its sweeps |
| sweep jobs | per unit {ward, average} × {nearest_centroid, wann, tacco, dot, tangram} → 1930 sweeps |
| settings | full K sweep, Leiden resolution 3.0, segmentation mask where the slice has one |
| hardware | TinyFat `work`, 8 cores, 64 GB, all thread pools pinned to 8; Tangram on CPU |

## Measurement

- `run_one.py` samples the RSS of the whole process tree (aligner children included) every 0.25 s.
- spatialaim's `timing.csv` times the scaffold and, per K, `map` and `analysis`.
- Joined into `phases.csv` (seconds + peak memory per phase and K) and `resources.json`
  (`wall_s` = `setup_s` + `sweep_s` + `rest_s`, `cache_hit`, peaks, host, spatialaim version).

## Files

| File | Role |
|---|---|
| `00_setup_server.sh` | envs: `spatial-aim==$SPATIALAIM_VERSION` from PyPI, tangram/tacco/dot envs, DOTr; `check` runs the checks |
| `01_scaffold.tinyfat.sbatch` | scaffold job, array index = the pair's Leiden unit |
| `02_sweep.tinyfat.sbatch` | sweep job, array index = unit, argument = mapper |
| `02b_sweep_by_k.tinyfat.sbatch` | one K of a sweep too long for 24 h, array index = K |
| `03_status.py` | done / failed / missing, cache misses, resubmit lines |
| `04_summarize.py` | → `summary/` |
| `config.sh`, `run_job.sh` | shared settings; body of every job |
| `grid.py`, `units.tsv` | the units, run folders, walltimes, `--array` lists |
| `run_one.py` | the measured runner: `scaffold`, `sweep`, `k`, `merge`, `levels` |

Every run is idempotent: `DONE` = skip, `FAILED` or no marker = rerun on resubmit.

## On the cluster

```
$WORK/spatialaim_exp1/
├── code/experiment1/           this folder; submit from here (logs/ = Slurm output)
├── data/01_Datasets_Paper/     incl. ST/Segmentation/; never move it during a run
├── envs/  conda_pkgs/  cache/
└── test/  run/                 run roots:
    ├── scaffolds/<pair_dir>/
    ├── results/<pair_dir>/<start>/<tree>/<mapper>/   k_comparison.csv, k_selection.json, k_levels.zip,
    │                                                 phases.csv, memory.tsv.gz, resources.json, run.log, DONE
    ├── split_k/                only for sweeps split by K
    ├── logs/ledger.tsv
    └── summary/
```

## Steps

Local commands are PowerShell, cluster commands bash on tinyx. SSH hosts `csnhr` (transfer) and
`tinyx` (ProxyJump `csnhr`) must be in `~/.ssh/config`.

**0. Release.** The spatial-aim version to pin must be on PyPI (publish a GitHub Release of
SpatialAIM-Repository).

**1. Clean slate and upload** (tinyx, then local):

```bash
rm -rf $WORK/spatialaim_exp1
```

```powershell
ssh tinyx "mkdir -p <WORK>/spatialaim_exp1/data <WORK>/spatialaim_exp1/code"
scp -r C:\Users\zi69hebi\Dev\10_Alignment\Data\01_Datasets_Paper csnhr:<WORK>/spatialaim_exp1/data/
scp -r C:\Users\zi69hebi\Dev\10_Alignment\SpatialAIM-Experiments\experiment1 csnhr:<WORK>/spatialaim_exp1/code/
```

**2. Envs** (tinyx). All checks must say `OK`, `spatialaim validate` no errors.

```bash
cd $WORK/spatialaim_exp1/code/experiment1
SPATIALAIM_VERSION=<version> bash 00_setup_server.sh
```

**3. One-pair test** (units 0 and 1 of pair 1: 1 scaffold + 20 sweeps):

```bash
source config.sh && conda activate spatialaim_env
export EXP1_ROOT=$EXP1_BASE/test
sbatch.tinyfat --array=0 01_scaffold.tinyfat.sbatch test
# once the scaffold is done:
for m in nearest_centroid wann tacco tangram; do
  sbatch.tinyfat -J exp1_sweep_$m --time=02:00:00 --array=0-1 02_sweep.tinyfat.sbatch test $m
done
sbatch.tinyfat -J exp1_sweep_dot --time=08:00:00 --array=0 02_sweep.tinyfat.sbatch test dot
sbatch.tinyfat -J exp1_sweep_dot --time=04:00:00 --array=1 02_sweep.tinyfat.sbatch test dot
python 03_status.py --failed                # until 21 / 21 done
python 04_summarize.py
```

Check: 0 failed, 0 cache misses; `meta.json` shows only the pinned version; `best_per_run.csv` has
30 rows; aligner rows have `peak_rss_children_mb` > 0; `resources_per_k.csv` is filled; no
tracebacks in `logs/`.

**4. Full grid.** Submit the sweeps only once all 94 scaffolds are done.

```bash
export EXP1_ROOT=$EXP1_BASE/run
sbatch.tinyfat --array=$(python grid.py array leiden)%20 01_scaffold.tinyfat.sbatch run
# once 03_status.py shows 94 scaffolds done:
for m in nearest_centroid wann tacco tangram; do
  sbatch.tinyfat -J exp1_sweep_$m --time=02:00:00 --array=0-192%20 02_sweep.tinyfat.sbatch run $m
done
sbatch.tinyfat -J exp1_sweep_dot --time=04:00:00 --array=$(python grid.py array annotation)%20 02_sweep.tinyfat.sbatch run dot
sbatch.tinyfat -J exp1_sweep_dot --time=08:00:00 --array=$(python grid.py array leiden)%20 02_sweep.tinyfat.sbatch run dot
```

**5. Leftovers.** `python 03_status.py --failed` prints one resubmit line per array; add options
before the script name, e.g. `--mem=128G` after an OOM or `--time=24:00:00`.

**5b. Sweeps over 24 h** (so far DOT on unit 189, breast cancer, Leiden), per tree without `DONE`:

```bash
sbatch.tinyfat --array=$(python run_one.py levels 189)%20 02b_sweep_by_k.tinyfat.sbatch run 189 dot ward
python run_one.py merge 189 dot ward        # once all K are done
```

Its times are sums over the K jobs, so `setup_s` is paid once per K; per-K `map`/`analysis` are unaffected.

**6. Summary and pull:**

```bash
python 04_summarize.py                      # expects missing: 0
```

```powershell
scp -r "csnhr:<WORK>/spatialaim_exp1/run/summary" C:\Users\zi69hebi\Dev\10_Alignment\Data\02_Results\experiment1
```

**7. Figures** (local): `python sections/02_reconstruction/panels.py`, then
`python sections/04_biology/panels.py --mapper <m>` per mapper; `make_excalidraw.py` only to
regenerate a drawing (it overwrites hand edits).

## `summary/`

| File | Content |
|---|---|
| `best_per_run.csv` | per sweep the compared K: `spatialaim` (Leiden, best overall), `baseline` (annotation, K = T), `annotation_best` (annotation, best overall), with the plotted metrics and `shape` |
| `all_k_per_run.csv.gz` | every `k_comparison.csv` row with the run's keys |
| `resources_sweep.csv` | per sweep: time split, peaks, `cache_hit`, plus `full_pipeline_s` (with the scaffold build) and `one_annotation_s` (annotation at K = T only) |
| `resources_per_k.csv` | per sweep and K: seconds and peak memory of `map` and `analysis` |
| `resources_scaffold.csv` | per pair: the scaffold build (`scaffold_s`, `scaffold_mb`) |
| `missing.csv`, `meta.json` | runs without `DONE`; counts and spatialaim version(s) |
