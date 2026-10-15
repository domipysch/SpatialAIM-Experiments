# Section 04 — Biological and topological plausibility (Fig. 2.4)

All data comes from `experiment1/`, reduced to one tidy table by `spatialaim_experiments.exp1`.

| Row | Criterion | Measures |
|---|---|---|
| a | spatial organisation | neighbourhood enrichment (self-diagonal z), local spatial purity z |
| b | transcriptional coherence | modularity over the slice's spots, modularity over the reference's cells |

Chance is 0 for all four by construction: the spatial pair are z-scores against a
spot-label shuffle, and modularity is the excess over a degree-preserving random
partition. 

SpatialAIM is, per pair, tree and mapper, the following configuration:
- take from the Leiden start the mapping at K = `best.overall`
- take from the annotation start(s) mapping at K = `best.overall`
- out of those, select the mapping with the best overall score

The annotation it is compared with is at K = T. The cost panels(d, e) time the Leiden runs, the full pipeline.

The nearest-centroid version opens with a top row the others do not have: the
win counts (one bar per metric) next to a worked example, the kidney mapped as
annotated and by SpatialAIM — see `04_kidney_example/README.md`.

## Running it

One plot per figure via `--mapper`.

```bash
python sections/04_biology/panels.py --mapper wann             # panels + panels.json
python sections/04_biology/make_excalidraw.py --mapper wann    # create excalidraw figure
python sections/04_biology/panels.py --all-mappers              # the same, aggregated over all five
python sections/04_biology/make_excalidraw.py --all-mappers
```

For nearest centroid, run the kidney example first (`04_kidney_example/01_run_mappings.py`,
then `02_panels.py`); `make_excalidraw.py` stops with a hint if its maps are missing.
