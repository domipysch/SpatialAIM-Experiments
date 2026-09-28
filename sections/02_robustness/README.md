# Section 02 — SpatialAIM is independent of the annotation's granularity (Fig. 2.2)

## Running it

Three steps: build the perturbed reference, run every mapper on the original and the perturbed
reference, draw the figure from those runs. Inputs are `Data/01_Datasets_Paper`; run from this
folder in PowerShell.

```powershell
$E = "C:\Users\zi69hebi\AppData\Local\miniforge3\envs\spatialaim_env"
$env:PATH = "$E;$E\Library\mingw-w64\bin;$E\Library\usr\bin;$E\Library\bin;$E\Scripts;$E\bin;$env:PATH"

# 1. the perturbed reference (once; seeded, so a re-run gives the same file)
#    -> Data/01_Datasets_Paper/Robustness/{split.h5ad, manifest.json}
& "$E\python.exe" 00_perturb.py

# 2. the sweeps: 3 runs per mapper (below), all five mappers; skips every run that
#    already has a k_selection.json, so an interrupted call just continues
powershell -File 01_run_robustness_experiments.ps1 -All

# 3. the figure, per mapper (analysis -> panels -> drawing -> PDF)
powershell -File 02_make_figures.ps1 -All
```

The three runs of step 2, all average linkage:

| run | reference | start clusters | K |
|---|---|---|---|
| `original_annotation` | `01_mouse-ssp` | annotation `cellType` (3 types); the unperturbed side of the invisible split and the annotation side of the granularity half | 3..1 |
| `invisible_split` | `Robustness/split.h5ad` | annotation `cellType_split` (4 types) | 4..1 |
| `original_leiden_overcluster` | `01_mouse-ssp` | Leiden over-clustering | L..1 |

Results land in `Data/02_Results/robustness/<run>/<mapper>/`, 
the analysis numbers in `Data/02_Results/robustness/_analysis/<mapper>/`. 

## The scripts

| File | What it does |
|------|--------------|
| `00_perturb.py` | halves Glutamatergic into `_A`/`_B` by nearest neighbour in panel space, then rewrites ten silent non-panel genes in the B half only. The 33 shared genes are never touched, so the split is invisible to the slice. → `Robustness/split.h5ad` + `manifest.json` (half sizes, panel centroid cosine, gene pairs) |
| `01_run_robustness_experiments.ps1` | the three sweeps for one mapper (`-Mapper`) or all (`-All`); skips finished runs, one log per run in `_logs/` |
| `02_make_figures.ps1` | runs 02a → 02b → 02c, per mapper (`-Mapper <name>` or `-All`) |
| `02a_analyse_invisible_split.py` | from the `invisible_split` run → `robustness/_analysis/<mapper>/results.json`: the spot share of each phantom half at K=4 (quoted in the text), plus SpatialAIM's proposed K and the score curves that `02b_panels.py` draws |
| `02b_panels.py` | the nine panels of the assembled figure + `panels.json` → `Data/02_Results/robustness/plots/<mapper>/`, copied into the vault folder the drawing embeds (`Excalidraw/figure_2/<mapper>/panels/`); also holds the run access (spot → type/state per K), the palette and the fine-annotation reading of the over-clustering states |
| `02c_make_excalidraw.py` | the drawing: letters, headings, the four arrows, three legends |

