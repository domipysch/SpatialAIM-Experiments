# Figure 4 — the kidney example

Figure 4 argues with numbers that SpatialAIM's mappings are more plausible than
mapping the author's annotation. This folder adds the picture: one tissue mapped
both ways, side by side, so a reader sees what a higher spatial assortativity
looks like.

The example is pair 4 (`04_mouse-kidney`, one MERFISH section of 126,547 cells)
with nearest centroid and average linkage. It was chosen as the reference where
SpatialAIM's gain in neighbourhood enrichment was largest, the spatial metric before
the assortativity. **Re-check that choice, and the numbers below, against the
Experiment 1 rerun**: the old values (z-scores) no longer apply.

## What is shown

Two mappings of Experiment 1 for this pair, mapper and tree, read from its summary
(`spatialaim_experiments.exp1`):

| side | start | K | kidney |
|---|---|---|---|
| annotation | by default the annotation `exp1.paired` compares against (of the reference's annotations, the one with the highest `score_overall`); `--ann_key` picks another | T, its number of types (identity cut) | **`cell_ontology_relabel`**, Liu et al.'s 13 types, K = 13 (`--ann_key cell_ontology_relabel`) |
| SpatialAIM | the start `exp1.best_start` picks | its `best.overall` | the Leiden over-clustering, K = 6 |

The kidney reference carries two annotations: TMS's own `cell_ontology_class` (22
classes) and Liu et al.'s coarse-grained `cell_ontology_relabel` (13). The figure shows
Liu's 13 types, the annotation of the paper the MERFISH section comes from. This
departs from the scatters of figure 4, whose kidney baseline is the 22 TMS classes:
they scored higher, by a small margin, in the previous Experiment 1 round
(`score_overall` 0.420 vs 0.419).

The per-spot mappings of Experiment 1 stay on the cluster, so step 1 recomputes the
two mappings locally, at those two Ks only, and checks their metrics against the
summary.

The figure is about spatial assortativity, not about which cell types these are,
so no state is named:

- **Maps** (whole section, small): every spot in the colour of its state. They only
  locate the zoom window, which is boxed.
- **Zooms** (large): the window with the neighbourhood graph drawn. An edge joining
  two spots of the same state is drawn in that state's colour, an edge between two
  states pale grey underneath. Assortativity is exactly the share of same-state
  edges, corrected for what randomly placed labels would give, so more colour and
  less grey is a higher assortativity.
- **Beside each map**: the global assortativity r on the window. The whole-section
  values and the per-cell-type means stay in `panels.json`.

The graph is SpatialAIM's own: squidpy's spatial kNN with 6 neighbours, the graph
the analysis step builds (`K_SPATIAL`). The window's assortativity uses SpatialAIM's
own function (`spatialaim.metrics.topology.assortativity`) on the graph rebuilt from
the window's spots. Unlike the old z-scores it does not grow with the number of
spots, so the window and the whole section are on one scale.

State colours (edges only): every annotated type has one colour. A SpatialAIM state
takes the colour of the annotated type that dominates its reference cells; a second
state dominated by the same type gets that colour's hue turned at the same luminance.

The window is set per pair in `ZOOM_WINDOWS` (`02_panels.py`). The kidney's, 9,000 x
7,875 units with about 4,000 spots, sits in the lower part of the section and spans
medulla, the NAPSA-high proximal-tubule ring and cortex (about 22 / 48 / 23 % of
its spots in SpatialAIM's states).

## Running it

From `spatialaim_env`, in this order:

```bash
python sections/04_biology/04_kidney_example/01_run_mappings.py --ann_key cell_ontology_relabel   # the two mappings -> Data/02_Results/biology_example/
python sections/04_biology/04_kidney_example/02_panels.py         # maps and zooms -> figure_4/nearest_centroid/kidney_example/
python sections/04_biology/panels.py --mapper nearest_centroid
python sections/04_biology/make_excalidraw.py --mapper nearest_centroid
```

Both scripts take `--pair`, `--mapper` and `--tree` (defaults: 4, `nearest_centroid`,
`average`). Step 1 skips a mapping whose K folder already exists; `--force` recomputes it.
It checks only the K, not the data, so pass `--force` whenever the dataset or the
Experiment 1 summary has changed.

`02_panels.py` needs squidpy for the graph: as for `spatialaim run`, prepend the
env's folders to PATH first, or it crashes silently. Its log names each SpatialAIM
state by its dominant annotated types (`SHORT_NAMES`), for reference only.

## Files

| file | does |
|---|---|
| `01_run_mappings.py` | decides the two mappings from the Experiment 1 summary, runs `spatialaim run` for each at its one K, writes `selection.json`, compares the metrics |
| `02_panels.py` | writes the two small maps, the two zooms with their graph and `panels.json` (sizes, metrics on the section and the window, zoom window) for `make_excalidraw.py` |
