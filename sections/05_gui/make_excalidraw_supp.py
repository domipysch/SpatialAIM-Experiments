# Assemble the supplementary figures of the GUI use case (Results 2.5) as
# Excalidraw drawings, with figure 5's primitives (make_excalidraw.py):
#
#   S10  the macrophage substates: the paper's markers (dotplot) and the Visium HD
#       check (registration, SPP1+ / SELENOP+ bins over cell types 2 / 3; the
#       per-bin counts and p-values go to the caption) -- panels from panels_s10.py
#       (figure_S10/panels/)
#   S7  choosing K in the GUI: the Choose K card and its trade-offs modal
#   S8  inspecting a result: report cards, cell-type modal (shape features),
#       confidence plot, the reference page's marker genes
#       (S7 / S8: screenshots_s7.py, figure_S7|S8/screenshots/)
#   S9 the criterion ablation over Experiment 1's Leiden sweeps: the K each
#       criterion set selects, each criterion's score range over K -- panels from
#       ablation_criteria/02_panels.py (figure_S9/panels/)
#
#   S11 neighbourhood enrichment of the Level1 annotation mapped directly: T cells,
#       endothelial cells and fibroblasts are mutual neighbours -- panel from
#       panels_s11.py (figure_S11/panels/)
#
#   python make_excalidraw_supp.py   -> figure_S7|S8|S9|S10|S11.excalidraw.md

import json
import logging
from string import ascii_lowercase

from spatialaim_experiments import paths
from spatialaim_experiments.excalidraw import WIDTH, Drawing

import make_excalidraw as fx
from make_excalidraw import GAP, PAD, SPAN, aspect, label, next_top, shot

logger = logging.getLogger(__name__)

S10 = paths.panel_dir("S10")
S7 = paths.figure_dir("S7") / "screenshots"
S8 = paths.figure_dir("S8") / "screenshots"
S9 = paths.panel_dir("S9")
S11 = paths.panel_dir("S11")


def row(m: dict, top: float, items, folder, frame: bool, h: float | None = None) -> float:
    """``items`` = (name, heading) side by side at one height: ``h``, or the one
    that fills SPAN."""
    if h is None:
        h = (SPAN - GAP * (len(items) - 1)) / sum(aspect(n, folder) for n, _ in items)
    x = PAD
    for name, heading in items:
        label(heading, x, top)
        w, _ = shot(m, name, x, top, height=h, folder=folder, frame=frame)
        x += w + GAP
    return h


def figure_s10(m: dict) -> float:
    top = PAD + fx.LETTER_DROP
    # the dotplot alone, at about 0.57 of the width (its print size, not stretched)
    h = row(m, top, [("dotplot", "Oliveira et al.'s macrophage markers, cell states 2 and 3")],
            S10, False, h=0.57 * SPAN / aspect("dotplot", S10))
    top = next_top(top + h)
    h = row(m, top, [("registration", "Visium HD registered onto Xenium"),
                     ("bins_spp1", "SPP1+ bins and cell state 2"),
                     ("bins_selenop", "SELENOP+ bins and cell state 3")], S10, False)
    return top + h + PAD


def rows(m: dict, folder, layout) -> float:
    top = PAD + fx.LETTER_DROP
    for i, items in enumerate(layout):
        h = row(m, top, items, folder, False)   # each screenshot draws its card outline
        top = next_top(top + h) if i < len(layout) - 1 else top + h
    return top + PAD


def figure_s7(m: dict) -> float:
    return rows(m, S7, ([("choosek", "Choose K card: Criterion scores across K")],
                        [("tradeoffs", "Criterion trade-offs (Finding Pareto-optimal K)")]))


def figure_s8(m: dict) -> float:
    return rows(m, S8, ([("report", "Report cards"),
                         ("modal_shape", "Inspecting cell states: cell shape features")],
                        [("confidence", "Mapping confidence"),
                         ("ref_markers", "Single-cell reference: Finding marker genes")]))


def figure_s9(m: dict) -> float:
    top = PAD + fx.LETTER_DROP
    h = row(m, top, [("selected_k_per_reference", "K selected by each criterion set")], S9, False)
    top = next_top(top + h)
    h = row(m, top, [("score_range_per_reference", "Range of each criterion score over K")],
            S9, False)
    return top + h + PAD


def figure_s11(m: dict) -> float:
    """The one heatmap, no heading (the caption names it), at its print width."""
    top = PAD
    _, h = fx.shot(m, "nhood_level1", PAD, top, width=0.5 * SPAN, folder=S11, frame=False)
    return top + h + PAD


def build(name: str, panels_dir, layout) -> None:
    fx.drawing = Drawing(f"f{name.lower()}", panels_dir, paths.VAULT)
    fx.letters = iter(ascii_lowercase)
    manifest_path = panels_dir / "panels.json"
    m = (json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.exists()
         else {"scale": 300, "panels": {}})
    height = layout(m)
    fx.drawing.check()
    out = paths.drawing(name)
    fx.drawing.write(out)
    logger.info("%s: canvas %.0f x %.0f units -> %s", name, WIDTH, height, out.name)


def main() -> None:
    build("S7", S7, figure_s7)
    build("S8", S8, figure_s8)
    build("S9", S9, figure_s9)
    build("S10", S10, figure_s10)
    build("S11", S11, figure_s11)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
