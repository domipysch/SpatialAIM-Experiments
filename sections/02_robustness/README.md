# Section 02 — AIM is independent of the annotation's granularity (Fig. 2.2)

Two experiments on **pair 0** (scRNA `01_mouse-ssp` × osmFISH `01_mouse-ssp`,
4839 cells, 33 genes, all shared), both with average linkage:

| | start clusters | what it asks |
|---|---|---|
| invisible split | the annotation, perturbed | an annotation **finer** than the panel can resolve — does AIM merge the phantom split back? |
| granularity | the default Leiden over-clustering | an annotation **coarser** than the panel supports — what does AIM propose instead? |

The slice is never touched, so the correct behaviour is known in advance.

## Running it

```powershell
# 1. build the perturbed reference (once) -> Data/01_Datasets/Robustness/split.h5ad
python perturb.py

# 2. the sweeps, per mapper; skips whatever already has a k_selection.json
powershell -File run_robustness.ps1 -Mapper tangram
powershell -File run_robustness.ps1 -Mapper dot -Only granularity

# 3. the figure, per mapper (analysis -> panels -> drawing -> PDF)
powershell -File make_figures.ps1 -Mapper tangram
powershell -File make_figures.ps1 -All
```

Results land in `Data/02_Results/robustness/{invisible_split,granularity}/<run>/<mapper>/`.
Reference aligners cost one `conda run` per K: the invisible-split runs are 3–4 K
levels (about a minute), the 44-level over-clustering sweep is ~18 minutes for DOT.

## The scripts

| File | What it does |
|------|--------------|
| `perturb.py` | halves Glutamatergic into `_A`/`_B` by nearest neighbour in panel space, then rewrites ten silent non-panel genes in the B half only. The 33 shared genes are never touched, so the split is invisible to the slice. |
| `run_robustness.ps1` | both experiments for one mapper, resumable, ledger in `_logs/` |
| `analyse_invisible_split.py` | the numbers the text quotes → `invisible_split/_analysis/<mapper>/{results.json,metrics.csv}`: how invisible the split is, what the mapper does with it at K=4, and whether K=3 reproduces the unperturbed map |
| `plot_invisible_split.py` | the first half's data access (run directories, spot→type, palette) **and** a standalone version of that half's figure |
| `plot_granularity.py` | the same for the second half: state composition against the fine annotation, the criterion curves |
| `panels.py` | the nine panels of the assembled figure + `panels.json` |
| `make_excalidraw.py` | the drawing: letters, headings, the four arrows, three legends |

The two `plot_*.py` scripts are imported by `panels.py` for their loaders and
colours; their own `main()` draws the older standalone figure that the assembled
one superseded. It is still a one-command look at the same data.

## Per mapper

Every script takes `--mapper` (default `tangram`, the figure the paper prints).
Any other mapper writes to `figure_2_2_<mapper>/` and
`figure_2_2_<mapper>.excalidraw.md`, so a copy can never overwrite the figure.

What the copies showed (September 2026):

* every mapper re-merges the invisible split at K=3 — **except DOT**, which
  proposes K=4 on every criterion, i.e. keeps the phantom split. Panel (g)'s
  heading follows the proposal, so the DOT copy reads "Split halves kept".
* "the division carries no spatial structure" is Tangram-specific. Nearest
  centroid divides the tissue into visible bands (local purity z = 48), DOT gives
  one half 71 % of the spots (z = 24), WANN and TACCO are mildly structured.
