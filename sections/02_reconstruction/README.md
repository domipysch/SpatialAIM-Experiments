# Section 02 — Reconstruction and granularity against the annotation (Fig. 2.2)

All data comes from `experiment1/`, reduced to one tidy table by `spatialaim_experiments.exp1`.

## Additional information

SpatialAIM is, per pair, tree and mapper, the following configuration:
- take from the Leiden start the mapping at K = `best.overall`
- take from the annotation start(s) mapping at K = `best.overall`
- out of those, select the mapping with the best overall score

The annotation it is compared with is at K = T. The cost panels(d, e) time the Leiden runs, the full pipeline.

Absolute cosine similarity is not comparable across references, so every panel is
a paired comparison *within* a reference, aggregated per reference.

## Running it

```bash
python sections/02_reconstruction/panels.py          # four panels + panels.json
python sections/02_reconstruction/make_excalidraw.py # create excalidraw
```

It reads `Data/02_Results/experiment1/summary/best_per_run.csv`, pulled from the cluster
(`experiment1/README.md`, step 6). `SPATIALAIM_EXP1` points at another pulled
result, e.g. the one-pair test.

## The panels

| Panel | |
|---|---|
| `delta_spot`, `delta_gene` | Δ cosine similarity, SpatialAIM − annotation, on a scale whose near-white band is stretched to ±0.1: almost every cell lies inside it, which is the row's point |
| `delta_k` | Δ number of cell types it took, clipped at ±20 (MTG runs to +69) |
| `mappers_cell`, `mappers_gene` | panels b, c: the five mappers at their own proposal, centred within each reference — the only way five mappers compare across ten incommensurable references; per cell and per gene, on one shared y axis labelled once in b |
| `umaps/{start,ward,average}.png` | row f–h: screenshots taken by hand from the GUI, not written by `panels.py`; kept in `Excalidraw/figure_2/umaps/` and placed by `make_excalidraw.py` at the row width |

Rows of the heat maps are mapper × linkage; their columns are the ten references.
The mapper names are drawn in Excalidraw, not in the images, so the three heat
maps keep identical geometry and can be read against each other.
