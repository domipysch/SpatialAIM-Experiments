# Section 04 — Biological and topological plausibility (Fig. 2.4)

Same corpus and same comparison as section 03 — `experiment1/` through
`aim_experiments.exp1`, AIM against the annotation baseline with the same mapper
on both sides — but on the other two criteria:

| Row | Criterion | Measures |
|---|---|---|
| a | spatial organisation | neighbourhood enrichment (self-diagonal z), local spatial purity z |
| b | transcriptional coherence | modularity over the slice's spots, modularity over the reference's cells |

Chance is 0 for all four by construction: the spatial pair are z-scores against a
spot-label shuffle, and modularity is the excess over a degree-preserving random
partition. Every run is above it, which the caption states.

Each row opens with a **sketch** of what its criterion means, drawn on a toy
graph with the same labels twice — observed above, permuted below. That is
exactly the null the numbers are measured against, so the sketch defines the
criterion and shows its chance level in one picture.

## Running it

```bash
python sections/04_biology/explore.py                # the numbers, printed
python sections/04_biology/panels.py                 # seven panels + panels.json
python sections/04_biology/make_excalidraw.py        # OVERWRITES the drawing
python -m aim_experiments.render 2_4                 # -> AIM-Overleaf/figures/fig_2_4.pdf
```

## One mapper per figure

The scatters show **one** mapper, `MAPPER` in `panels.py` (`nearest_centroid`).
Both axes of a panel run the same aligner, so the mapper is a paired factor and
never a confound; all five at once turned each panel into 50 points for 10
datasets and read as clutter. That the result holds for every aligner belongs in
the text.

`--mapper <other>` writes a full copy of the figure into `figure_2_4_<mapper>/`
and never touches the figure's own panels. `--all-mappers` is a trial: 50 points
per panel, preview folder only, no manifest.

Above the parity diagonal, out of ten references (September 2026):

| mapper | nhood enrichment | local purity | modularity spatial | modularity reference |
|---|---|---|---|---|
| nearest centroid | 7 | 8 | 8 | 4 |
| Tangram | 7 | 5 | 6 | 4 |
| WANN | 7 | 5 | 6 | 3 |
| TACCO | 7 | 6 | 7 | 5 |
| DOT | 6 | 4 | 8 | 6 |

Modularity over the reference's cells is the one panel where AIM usually is not
ahead, and it is also the one whose Δ correlates positively with Δ number of cell
types (ρ ≈ +0.4 to +0.8, negative for the other three) — it mostly reflects AIM
using fewer cell types than the annotation.
