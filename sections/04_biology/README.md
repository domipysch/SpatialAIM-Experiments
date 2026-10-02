# Section 04 — Biological and topological plausibility (Fig. 2.4)

All data comes from `experiment1/`, reduced to one tidy table by `spatialaim_experiments.exp1`.

| Row | Criterion | Measures |
|---|---|---|
| a | spatial organisation | neighbourhood enrichment (self-diagonal z), local spatial purity z |
| b | transcriptional coherence | modularity over the slice's spots, modularity over the reference's cells |

Chance is 0 for all four by construction: the spatial pair are z-scores against a
spot-label shuffle, and modularity is the excess over a degree-preserving random
partition. 

## Running it

One plot per figure via `--mapper`.

```bash
python sections/04_biology/panels.py --mapper wann             # seven panels + panels.json
python sections/04_biology/make_excalidraw.py --mapper wann    # create excalidraw figure
python sections/04_biology/panels.py --all-mappers              # the same, aggregated over all five
python sections/04_biology/make_excalidraw.py --all-mappers
```
