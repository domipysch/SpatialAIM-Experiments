# Section 03 — Reconstruction and granularity against the annotation (Fig. 2.3)

No runs of its own: everything comes from `experiment1/`, reduced to one tidy
table by `spatialaim_experiments.exp1`.

The comparison, in every panel: **SpatialAIM** (Leiden over-clustering, K = `best.overall`)
against **the annotation** (the identity cut at K = T of the same reference),
joined per pair, linkage and mapper, so the same aligner runs on both sides and
only the starting partition differs. Where a reference carries two annotations,
the baseline is the one whose identity cut scores highest *overall* — the same
harmonic mean by which SpatialAIM picks its own K.

Absolute cosine similarity is not comparable across references (33-gene panels
against 981-gene ones, and provenance sets a per-pair ceiling), so every panel is
a paired comparison *within* a reference, aggregated per reference.

## Running it

```bash
python sections/03_reconstruction/panels.py          # four panels + panels.json
python sections/03_reconstruction/make_excalidraw.py # OVERWRITES the drawing
python -m spatialaim_experiments.render 2_3                 # -> SpatialAIM-Overleaf/figures/fig_2_3.pdf
```

The first call builds `outputs/experiment1/runs.csv` if it is missing — about a
minute over 2630 run roots. Delete that file (or `python -m spatialaim_experiments.exp1
--force`) after collecting new runs.

## The panels

| Panel | |
|---|---|
| `delta_spot`, `delta_gene` | Δ cosine similarity, SpatialAIM − annotation, on a scale whose near-white band is stretched to ±0.1: almost every cell lies inside it, which is the row's point |
| `delta_k` | Δ number of cell types it took, clipped at ±20 (MTG runs to +69) |
| `mappers` | the five mappers at their own proposal, centred within each reference — the only way five mappers compare across ten incommensurable references |

Rows of the heat maps are mapper × linkage; their columns are the ten references.
The mapper names are drawn in Excalidraw, not in the images, so the three heat
maps keep identical geometry and can be read against each other.
