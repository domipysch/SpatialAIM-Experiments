"""
The panel of supplementary figure S11: neighbourhood enrichment (z-scores, 6-NN graph)
of CRC P2's Level1 annotation mapped directly (nearest centroid, K = 9), T cells,
endothelial cells and fibroblasts boxed -- the types SpatialAIM merges at K = 10.

The heatmap and the squidpy computation are analysis/nhood_level1.py; this
reuses its saved z-scores (Data/03_GUI/Analysis/nhood_level1/nhood_z.csv; run that script first) and writes
the panel to Excalidraw/figure_S11/panels/.

    python sections/05_gui/panels_s11.py
"""

from __future__ import annotations

import importlib.util
import logging
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import pandas as pd  # noqa: E402
from spatialaim_experiments import paths  # noqa: E402
from spatialaim_experiments.style import PanelWriter  # noqa: E402

FIGURE = "S11"
ANALYSIS = Path(__file__).resolve().parent / "analysis"   # the scripts; their outputs: paths.GUI_ANALYSIS


def nhood_level1():
    """nhood_level1.py as a module (its heatmap, not its squidpy run)."""
    spec = importlib.util.spec_from_file_location("nhood_level1", ANALYSIS / "nhood_level1.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    nh = nhood_level1()
    z = pd.read_csv(nh.OUT / "nhood_z.csv", index_col=0)
    panels = PanelWriter(paths.panel_dir(FIGURE), preview=paths.outputs("05_gui"))
    panels.save(nh.heatmap(z), "nhood_level1")
    panels.write_manifest()


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
