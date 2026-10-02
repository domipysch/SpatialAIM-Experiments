# SpatialAIM-Experiments

## Layout

```
experiment1/            the prerequisite sweep: every pair x start mode x tree x
                        mapper, run under Slurm on NHR TinyFat (Tangram on CPU) with time
                        and memory measured per run. Sections 03 and 04 read its
                        pulled summary.
sections/               one folder per results section of the paper, each holding
                        the runs it needs and the scripts that draw its figure
  01_overview/          Fig. 2.1  what SpatialAIM does            (hand-drawn)
  02_reconstruction/    Fig. 2.2  reconstruction & granularity vs the annotation
  03_robustness/        Fig. 2.3  annotation granularity   (own runs, pair 1)
  04_biology/           Fig. 2.4  biological & topological plausibility
src/spatialaim_experiments/    the helpers every section shares (installed, see below)
```

Each section has its own README with the order its scripts run in and what each one writes.

## Working directories

Only `src/spatialaim_experiments/paths.py` spells out an absolute path; everything else
goes through it, and every location has an environment override:

| What | Default | Override |
|------|---------|----------|
| datasets + results | `../Data` | `SPATIALAIM_DATA` |
| Experiment 1 summary | `../Data/02_Results/experiment1` | `SPATIALAIM_EXP1` |
| Excalidraw drawings | `../../00_Obsidian` | `SPATIALAIM_VAULT` |
| previews and caches | `./outputs` | `SPATIALAIM_EXPERIMENTS_OUT` |

Setting `SPATIALAIM_VAULT` to a scratch folder is the safe way to try a figure out: the
whole pipeline then writes there and no drawing in the vault is touched.

## How a figure is made

The same four steps in every section, and they are deliberately separate — the
expensive part is the runs, and the drawing is edited by hand afterwards:

1. **runs** — `spatialaim run` produces metrics under `Data/02_Results`. Section 02 has
   its own small runs; sections 03 and 04 read `experiment1/`.
2. **`panels.py`** — draws each panel as its own PNG into the Excalidraw folder
   in the vault, plus `panels.json` recording each panel's size in print inches.
   A panel carries only what is anchored to data.
3. **`make_excalidraw.py`** — assembles those panels into
   `figure_<n>.excalidraw.md`: letters, headings, arrows, legends.
4. **PDF** — exported by hand from the drawing in Obsidian (Excalidraw's own export)
   into the paper's `figures/fig_2_<n>.pdf`.

## Setup

The experiments run in `spatialaim_env`, the same conda environment the tool uses, which
already has the scientific stack:

```bash
conda activate spatialaim_env
pip install -e ../SpatialAIM-Repository --no-deps   # the method, as `spatialaim`
pip install -e . --no-deps                          # these helpers
```

`--no-deps` on both: every dependency is already in the environment, and it stops
pip from trying to fetch `spatial-aim` from an index instead of using the local
checkout. Check it:

## Conventions

* **Aggregate per reference, never per pair.** 72 of the 88 pairs are retina slices.
* **SpatialAIM's level is `best.overall`** from `k_selection.json`; the annotation
  baseline is the identity cut at K = T of the same run tree, with the same
  mapper on both sides.
