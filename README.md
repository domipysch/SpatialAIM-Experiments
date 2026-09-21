# SpatialAIM-Experiments

Everything that *uses* AIM: the experiments behind the paper, the analyses of
their output, and the figures.

The method itself is not here. It is the **`spatial-aim`** package in
[SpatialAIM-Repository](../SpatialAIM-Repository), which this repository imports
as `aim` and never modifies. The dependency runs one way only:

```
SpatialAIM-Repository   the tool     ->  installed as `spatial-aim`
SpatialAIM-Experiments  this repo    ->  imports `aim`, calls the `aim` CLI
```

Nothing in the tool repository reads anything from here.

## Layout

```
experiment1/            the prerequisite sweep: every pair x start mode x tree x
                        mapper, run on the compute servers. Sections 03 and 04
                        read its collected metrics.
sections/               one folder per results section of the paper, each holding
                        the runs it needs and the scripts that draw its figure
  01_overview/          Fig. 2.1  what AIM does            (hand-drawn)
  02_robustness/        Fig. 2.2  annotation granularity   (own runs, pair 0)
  03_reconstruction/    Fig. 2.3  reconstruction & granularity vs the annotation
  04_biology/           Fig. 2.4  biological & topological plausibility
supplement/             the supplementary dataset tables
src/aim_experiments/    the helpers every section shares (installed, see below)
outputs/                previews, caches, trial renderings -- all regenerable,
                        all git-ignored
```

Each section has its own README with the order its scripts run in and what each
one writes.

## Where things live

Only `src/aim_experiments/paths.py` spells out an absolute path; everything else
goes through it, and every location has an environment override:

| What | Default | Override |
|------|---------|----------|
| datasets + results | `../Data` | `AIM_DATA` |
| Excalidraw drawings | `../../00_Obsidian` | `AIM_VAULT` |
| the paper | `../AIM-Overleaf` | `AIM_PAPER` |
| previews and caches | `./outputs` | `AIM_EXPERIMENTS_OUT` |

Setting `AIM_VAULT` to a scratch folder is the safe way to try a figure out: the
whole pipeline then writes there and no drawing in the vault is touched.

## How a figure is made

The same four steps in every section, and they are deliberately separate — the
expensive part is the runs, and the drawing is edited by hand afterwards:

1. **runs** — `aim run` produces metrics under `Data/02_Results`. Section 02 has
   its own small runs; sections 03 and 04 read `experiment1/`.
2. **`panels.py`** — draws each panel as its own PNG into the Excalidraw folder
   in the vault, plus `panels.json` recording each panel's size in print inches.
   A panel carries only what is anchored to data.
3. **`make_excalidraw.py`** — assembles those panels into
   `figure_<n>.excalidraw.md`: letters, headings, arrows, legends.
4. **render** — `python -m aim_experiments.render 2_3` turns the drawing into
   `AIM-Overleaf/figures/fig_2_3.pdf`. Pass `--out something.png` for a look that
   never touches the paper.

> **`make_excalidraw.py` overwrites the drawing, hand edits included.** Step 2 is
> always safe (it replaces images the drawing embeds by link); step 3 is not,
> once the figure has been adjusted in Obsidian. Every drawing in the vault is
> currently hand-edited — re-run step 3 only when you mean to start that figure's
> layout over.

## Setup

The experiments run in `aim_env`, the same conda environment the tool uses, which
already has the scientific stack:

```bash
conda activate aim_env
pip install -e ../SpatialAIM-Repository --no-deps   # the method, as `aim`
pip install -e . --no-deps                          # these helpers
```

`--no-deps` on both: every dependency is already in the environment, and it stops
pip from trying to fetch `spatial-aim` from an index instead of using the local
checkout. Check it:

```bash
python -c "import aim, aim_experiments; print(aim.__file__)"
```

On Windows the conda environment has to be on `PATH`, not just its interpreter —
otherwise the native libraries load from the wrong place and a sweep dies with
`0xC0000409`. The `.ps1` drivers do this themselves; for a bare shell:

```powershell
$E = "C:\Users\zi69hebi\AppData\Local\miniforge3\envs\aim_env"
$env:PATH = "$E;$E\Library\mingw-w64\bin;$E\Library\usr\bin;$E\Library\bin;$E\Scripts;$E\bin;$env:PATH"
```

The reference aligners (Tangram, TACCO, DOT) each need their own conda env —
`tangram_env`, `tacco_env`, `dot_env` — created from the specs shipped inside the
`spatial-aim` wheel. See the tool repository's README, or `experiment1/README.md`
for the server-side version.

## Conventions

* **Aggregate per reference, never per pair.** 71 of the 92 pairs are retina
  islands; pooling over pairs is the retina's answer, not the corpus's.
* **The 92 whole pairs only** (`pairs_only_raw_counts_11_12_whole.csv`), not the
  244-pair database and not the tiled 11/12 variants.
* **AIM's level is `best.overall`** from `k_selection.json`; the annotation
  baseline is the identity cut at K = T of the same run tree, with the same
  mapper on both sides.
* Nothing here writes to `Data/01_Datasets` except `sections/02_robustness/perturb.py`,
  which builds one perturbed reference of its own.
