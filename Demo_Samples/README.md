# Demo results — on NHR@FAU TinyFat

Computes the online demo's precomputed results: every pair of the demo's `pairs.csv` × start
(Leiden + every `CellTypeKey*`) × linkage {average, ward} × mapper {nearest_centroid, wann,
tacco, tangram, dot}, full K sweep, Leiden resolution 3.0, segmentation mask where the slice has
one. Unlike Experiment 1 it keeps every file (the per-K mappings included) and writes them
straight into the layout the demo reads, so `results/` is downloaded as it is.

For the 10 demo pairs: 230 runs, ≈ 250 core-h × 8 in Experiment 1, ≈ 2 × 17 GB of results.

## Files

| File | Role |
|---|---|
| `00_setup_server.sh` | installs `spatial-aim==$SPATIALAIM_VERSION` into Experiment 1's `spatialaim_env`; `check` runs the checks |
| `01_scaffold.tinyfat.sbatch` | scaffold job, array index = PairID |
| `02_sweep.tinyfat.sbatch` | sweep job, array index = run (row of `runs.tsv`) |
| `02b_sweep_by_k.tinyfat.sbatch` | one K of a run too long for 2 h (`grid.SPLIT_BY_K`: every Experiment 1 run from 1.8 h, DOT and Tangram on 002, 003, 004, 093, 094), array index = K |
| `03_submit_sweeps.sh` | submits all sweeps, one 2 h array per class (fast, tangram, dot) |
| `config.sh`, `run_job.sh` | shared settings; body of every job |
| `grid.py`, `runs.tsv` | the runs, walltime classes, `--array` lists |
| `run_one.py` | `scaffold`, `sweep`, `k`, `merge`, `levels`, `status` |

Every run is idempotent: `DONE` = skip, `FAILED` or no marker = rerun on resubmit.

## On the cluster

```
$WORK/spatialaim_demo/
├── code/Demo_Samples/   this folder; submit from here (logs/ = Slurm output)
├── datasets/            = Data/04_Demo_sample/datasets (scRNA/, ST/, pairs.csv)
├── results/             the download: <linkage>/<pair_dir>/<start>/<mapper>/
└── state/               DONE/FAILED + run.log per run, split_k/; not downloaded
```

Envs are Experiment 1's in `$WORK/spatialaim_exp1/envs` (`DEMO_ENVS` in `config.sh`), used as
they are; the setup only upgrades spatial-aim in `spatialaim_env`.

## Steps

Local commands are PowerShell, cluster commands bash on tinyx. SSH hosts `csnhr` (transfer) and
`tinyx` (ProxyJump `csnhr`) as for Experiment 1.

**0. Release.** The spatial-aim version must be on PyPI and write what the deployed demo reads
(the same version as the demo image, or newer).

**1. Upload** (tinyx, then local):

```bash
mkdir -p $WORK/spatialaim_demo/code
```

```powershell
scp -r C:\Users\zi69hebi\Dev\10_Alignment\Data\04_Demo_sample\datasets csnhr:<WORK>/spatialaim_demo/
scp -r C:\Users\zi69hebi\Dev\10_Alignment\SpatialAIM-Experiments\Demo_Samples csnhr:<WORK>/spatialaim_demo/code/
```

**2. Env** (tinyx). The checks must say `OK`, `spatialaim validate` no errors.

```bash
cd $WORK/spatialaim_demo/code/Demo_Samples
SPATIALAIM_VERSION=<version> bash 00_setup_server.sh
```

**3. Scaffolds**, one per pair:

```bash
source config.sh && conda activate $SPATIALAIM_ENV
sbatch.tinyfat --array=$(python grid.py pairs) 01_scaffold.tinyfat.sbatch
python run_one.py status                     # until all scaffolds are done
```

**4. Sweeps:**

```bash
bash 03_submit_sweeps.sh
python run_one.py status                     # failed runs + one resubmit line per class
```

Add options before the script name when resubmitting, e.g. `--mem=128G` after an OOM or
`--time=24:00:00`. A run too long for its class goes into `grid.SPLIT_BY_K`; `status` prints its
lines.

**5. Merge the runs split by K**, once all their K are done:

```bash
python run_one.py merge                      # every split run; one with: merge <run>
```

**6. Download** once `status` shows every run done:

```powershell
scp -r "csnhr:<WORK>/spatialaim_demo/results" C:\Users\zi69hebi\Dev\10_Alignment\Data\04_Demo_sample\
```

Then copy it into the demo server's data folder (`/srv/spatialaim-demo/data/results/`). Pairs
show up there as soon as they have a finished start, so a partial download already works.
