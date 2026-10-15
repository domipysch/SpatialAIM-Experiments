# Figure 4 — the kidney example

Figure 4 argues with numbers that SpatialAIM's mappings are more plausible than
mapping the author's annotation. This folder adds the picture: one tissue mapped
both ways, side by side, so a reader sees what a higher neighbourhood enrichment
looks like.

The example is pair 4 (`04_mouse-kidney`, one MERFISH section of 126,547 cells)
with nearest centroid and average linkage — the reference where SpatialAIM's gain
is largest. Mapped as annotated (9 types), the medulla is a speckle of cortical
types (podocytes, distal convoluted tubule). SpatialAIM proposes 6 states, and the
section separates into cortex, a band at the cortex–medulla border and medulla.

## What is shown

Exactly the two mappings Experiment 1 compares for this pair, mapper and tree,
read from its summary (`spatialaim_experiments.exp1`) — nothing is chosen here:

| side | start | K |
|---|---|---|
| annotation | the annotation `exp1.paired` compares against (highest `score_overall`; for the kidney `celltype`) | T, its number of types (identity cut) |
| SpatialAIM | the start `exp1.best_start` picks (for the kidney the Leiden over-clustering) | its `best.overall` |

The per-spot mappings of Experiment 1 stay on the cluster, so step 1 recomputes the
two mappings locally, at those two Ks only, and checks their metrics against the
summary.

Colours: every annotated type has one colour. A SpatialAIM state takes the colour
of the annotated type that dominates its reference cells, so a merge shows as one
colour covering several types; a second state dominated by the same type gets
that colour's hue turned at the same luminance (and a numbered name in the key).
All colours are saturated and of similar weight: a pale colour in one map, or a
near-black in the other, would make one map look washed out against the other
for reasons of palette alone. Spots are drawn in a random order (fixed seed):
which colour ends up on top where spots overlap is then decided by how many
spots each state has there, not by its id.

Under each map, a zoom shows one window of the section (boxed on both maps) at
spot level, where a mixed and a coherent mapping differ visibly whatever the
colours. The window is set per pair in `ZOOM_WINDOWS` (`02_panels.py`); the
kidney's crosses the cortex–medulla border.

## Running it

From `spatialaim_env`, in this order:

```bash
python sections/04_biology/04_kidney_example/01_run_mappings.py   # the two mappings -> Data/02_Results/biology_example/
python sections/04_biology/04_kidney_example/02_panels.py         # the two maps -> figure_4/nearest_centroid/kidney_example/
python sections/04_biology/panels.py --mapper nearest_centroid
python sections/04_biology/make_excalidraw.py --mapper nearest_centroid
```

Both scripts take `--pair`, `--mapper` and `--tree` (defaults: 4, `nearest_centroid`,
`average`). Step 1 skips a mapping that already exists; `--force` recomputes it.

## Files

| file | does |
|---|---|
| `01_run_mappings.py` | decides the two mappings from the Experiment 1 summary, runs `spatialaim run` for each at its one K, writes `selection.json`, compares the metrics |
| `02_panels.py` | colours both mappings consistently and writes the two maps, the two zooms and `panels.json` (sizes, keys, metrics, zoom window) for `make_excalidraw.py` |
