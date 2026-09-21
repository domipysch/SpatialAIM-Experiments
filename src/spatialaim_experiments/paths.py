"""Every filesystem location the experiments touch, resolved in one place.

Nothing else in this repository spells out an absolute path. Each location has a
default that matches this workstation and an environment variable that overrides
it, so the same scripts run on another machine without editing a line:

    SPATIALAIM_DATA             the datasets and results tree      (Data/)
    SPATIALAIM_VAULT            the Obsidian vault holding the Excalidraw drawings
    SPATIALAIM_PAPER            the Overleaf project               (SpatialAIM-Overleaf/)
    SPATIALAIM_EXPERIMENTS_OUT  where previews and caches are written (outputs/)

The results this repository reads are produced by ``spatial-aim`` (the tool
repository) and live under ``SPATIALAIM_DATA``; nothing is ever written back there by a
figure script.
"""

from __future__ import annotations

import os
from pathlib import Path

# The repository root, found from this file: src/spatialaim_experiments/paths.py. It is
# correct for an editable install (the only one documented in the README) and is
# overridden by SPATIALAIM_EXPERIMENTS_OUT wherever that does not hold.
REPO = Path(__file__).resolve().parents[2]
WORKSPACE = Path(r"C:\Users\zi69hebi\Dev\10_Alignment")


def _env(name: str, default: Path) -> Path:
    value = os.environ.get(name)
    return Path(value).expanduser() if value else default


# ---------------------------------------------------------------------------
# inputs: the datasets and the results of the sweeps
# ---------------------------------------------------------------------------
DATA = _env("SPATIALAIM_DATA", WORKSPACE / "Data")
DATASETS = DATA / "01_Datasets"
RESULTS = DATA / "02_Results"

# The 92 pairs the paper uses -- the whole-slice variants of datasets 11 and 12,
# not the tiled ones. Every aggregation in this repository is restricted to it.
PAIRS_CSV = DATASETS / "pairs_only_raw_counts_11_12_whole.csv"
SC_INDEX = DATASETS / "scRNA" / "index.csv"
ST_INDEX = DATASETS / "ST" / "index.csv"

# the two result trees the figures read
EXPERIMENT1 = RESULTS / "experiment1"
ROBUSTNESS = RESULTS / "robustness"

# ---------------------------------------------------------------------------
# outputs: the drawings in the vault, the PDFs in the paper, previews here
# ---------------------------------------------------------------------------
VAULT = _env("SPATIALAIM_VAULT", Path(r"C:\Users\zi69hebi\Dev\00_Obsidian"))
EXCALIDRAW = VAULT / "Excalidraw"
# The Overleaf checkout was called AIM-Overleaf before the rename to SpatialAIM.
# Either name works, so a machine that has not renamed its folder yet still
# resolves; SPATIALAIM_PAPER overrides both.
_PAPER_NAMES = ("SpatialAIM-Overleaf", "AIM-Overleaf")
PAPER = _env("SPATIALAIM_PAPER",
             next((WORKSPACE / name for name in _PAPER_NAMES
                   if (WORKSPACE / name).exists()),
                  WORKSPACE / _PAPER_NAMES[0]))
PAPER_FIGURES = PAPER / "figures"
OUTPUTS = _env("SPATIALAIM_EXPERIMENTS_OUT", REPO / "outputs")


def outputs(section: str) -> Path:
    """This section's own folder for previews, caches and trial renderings.

    Everything here is regenerable and git-ignored; nothing the paper includes
    is written here except by an explicit ``--out``.
    """
    path = OUTPUTS / section
    path.mkdir(parents=True, exist_ok=True)
    return path


# ---------------------------------------------------------------------------
# one figure, by its number in the paper ("2_2"), optionally per mapper
# ---------------------------------------------------------------------------

def _stem(figure: str, mapper: str | None) -> str:
    return f"figure_{figure}" if mapper is None else f"figure_{figure}_{mapper}"


def panel_dir(figure: str, mapper: str | None = None) -> Path:
    """Where a figure's panel PNGs go, next to the drawing that embeds them."""
    return EXCALIDRAW / _stem(figure, mapper)


def drawing(figure: str, mapper: str | None = None) -> Path:
    """The Excalidraw document assembled from those panels."""
    return EXCALIDRAW / f"{_stem(figure, mapper)}.excalidraw.md"


def figure_pdf(figure: str, mapper: str | None = None) -> Path:
    """Where the paper includes the rendered figure from."""
    suffix = figure if mapper is None else f"{figure}_{mapper}"
    return PAPER_FIGURES / f"fig_{suffix}.pdf"
