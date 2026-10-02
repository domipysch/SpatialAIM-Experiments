"""Every filesystem location the experiments touch, resolved in one place.

Nothing else in this repository spells out an absolute path. Each location has a
default that matches this workstation and an environment variable that overrides
it, so the same scripts run on another machine without editing a line:

    SPATIALAIM_DATA             the datasets and results tree      (Data/)
    SPATIALAIM_EXP1             one Experiment 1 result, as pulled (Data/02_Results/experiment1)
    SPATIALAIM_VAULT            the Obsidian vault holding the Excalidraw drawings
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
DATASETS = DATA / "01_Datasets_Paper"
RESULTS = DATA / "02_Results"

# The paper's 88 pairs (10 references, 1-based PairIDs).
PAIRS_CSV = DATASETS / "pairs.csv"
SC_INDEX = DATASETS / "scRNA" / "index.csv"
ST_INDEX = DATASETS / "ST" / "index.csv"

# Experiment 1 as pulled from the cluster: its summary/ is all the figures read
# (experiment1/04_summarize.py). SPATIALAIM_EXP1 points elsewhere, e.g. at the
# one-pair test.
EXPERIMENT1 = _env("SPATIALAIM_EXP1", RESULTS / "experiment1")
EXP1_SUMMARY = EXPERIMENT1 / "summary"
ROBUSTNESS = RESULTS / "robustness"

# ---------------------------------------------------------------------------
# outputs: the drawings in the vault, previews here
# ---------------------------------------------------------------------------
VAULT = _env("SPATIALAIM_VAULT", Path(r"C:\Users\zi69hebi\Dev\00_Obsidian"))
EXCALIDRAW = VAULT / "Excalidraw"
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
# one results figure, by its number ("3" for Fig. 2.3), optionally per mapper;
# one folder per figure in the vault, a mapper's copy in a subfolder of it:
#
#   Excalidraw/figure_3/figure_3.excalidraw.md        the drawing
#   Excalidraw/figure_3/panels/                        the PNGs it embeds
#   Excalidraw/figure_3/dot/figure_3_dot.excalidraw.md the same figure for DOT
#   Excalidraw/figure_3/dot/panels/
# ---------------------------------------------------------------------------

def figure_dir(figure: str, mapper: str | None = None) -> Path:
    """A figure's folder in the vault, or its mapper copy's subfolder."""
    folder = EXCALIDRAW / f"figure_{figure}"
    return folder if mapper is None else folder / mapper


def panel_dir(figure: str, mapper: str | None = None) -> Path:
    """Where a figure's panel PNGs go, next to the drawing that embeds them."""
    return figure_dir(figure, mapper) / "panels"


def drawing(figure: str, mapper: str | None = None) -> Path:
    """The Excalidraw document assembled from those panels."""
    stem = f"figure_{figure}" if mapper is None else f"figure_{figure}_{mapper}"
    return figure_dir(figure, mapper) / f"{stem}.excalidraw.md"

