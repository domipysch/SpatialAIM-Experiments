"""
GUI screenshots of supplementary figures S7 and S8 (the GUI features figure 5 only
names), same run and setup as screenshots.py (GUI on 8503, headless Chrome on 9333),
Leiden start, nearest centroid, K = 10. 2 px per CSS px.

Excalidraw/figure_S7/screenshots/ -- choosing K:

    choosek.png       the open "Choose K" card: the criteria over K, best K per criterion
                      (cut above its divider)
    tradeoffs.png     its "Criterion trade-offs" modal (without the note below the plots)

Excalidraw/figure_S8/screenshots/ -- inspecting a result:

    report.png        the Report expander down to the mapping-confidence card (its
                      shape-features card is left out)
    modal_shape.png   cell type 2's modal, "Shape features" view
    confidence.png    the spatial plot with "Plot confidence" on, and its controls
    ref_markers.png   the Single-cell reference page's "Cell type marker genes" card
                      for the Level1 start at K = T = 9, i.e. the Level1 annotation
"""

from __future__ import annotations

import sys
import time

from PIL import Image
from spatialaim_experiments import paths

import screenshots as sc
from cards import as_card
from screenshots import DSF, W, click, p

OUT = {name: paths.figure_dir(fig) / "screenshots"
       for fig, names in (("S7", ("choosek", "tradeoffs")),
                          ("S8", ("report", "modal_shape", "confidence", "ref_markers")))
       for name in names}
# heights kept from the top (CSS px): the report down to the mapping-confidence
# card; the Choose K card down to its per-criterion best K (above the divider);
# the trade-offs modal down to its Export buttons (without the grey note)
REPORT_TOP, CHOOSEK_TOP, TRADEOFFS_TOP = 990, 505, 509
NO_TIPS = ("const s=document.createElement('style'); s.id='noTip';"
           "s.textContent='[role=tooltip],[data-baseweb=tooltip],[data-baseweb=popover]"
           "{display:none!important}'; document.head.appendChild(s)")


def expander(title: str) -> str:
    return (f"[...document.querySelectorAll('[data-testid=stExpander]')]"
            f".find(e=>{{const s=e.querySelector('summary'); return s && s.innerText.includes('{title}')}})")


def settle(scope: str, timeout: int = 180) -> None:
    """Wait until ``scope`` shows no spinner / "Computing" note, then a little more
    for Plotly to draw."""
    for _ in range(timeout // 3):
        busy = p.js(f"(()=>{{const e={scope}; return !e || e.innerText.includes('Computing') || "
                    f"e.innerText.includes('Ranking') || !!e.querySelector('[data-testid=stSpinner]')}})()")
        if not busy:
            break
        time.sleep(3)
    time.sleep(8)


def shot(js: str, name: str, height: int | None = None) -> None:
    p.js("document.activeElement && document.activeElement.blur();" + NO_TIPS)
    p.move(W - 5, 5); time.sleep(2)
    r = p.rect(js)
    h = r["h"] if height is None else min(r["h"], height)
    path = OUT[name.removesuffix(".png")] / name
    p.shot(str(path), clip=(r["x"], r["y"], r["w"], h))
    as_card(path)
    p.js("document.getElementById('noTip')?.remove()")
    print(f"wrote {name} ({r['w']:.0f} x {h:.0f} CSS px)")


def open_expander(title: str, height: int = 3600) -> str:
    p.viewport(W, height, DSF); time.sleep(4)
    exp = expander(title)
    click(f"{exp}.querySelector('summary')", wait=8)
    settle(exp)
    return exp


def main(only: set[str]) -> None:
    """``only``: the step names to take (choosek, report, modal, confidence,
    reference); empty = all."""
    for d in set(OUT.values()):
        d.mkdir(parents=True, exist_ok=True)
    sc.open_result()
    sc.collapse_sidebar()
    want = lambda step: not only or step in only  # noqa: E731

    if want("choosek"):
        choose_k()
    if want("report"):
        report()
    if want("modal"):
        modal_shape()
    if want("confidence"):
        confidence()
    if want("reference"):
        reference_markers()


def choose_k() -> None:
    """The Choose K card and, from it, the criterion trade-offs modal."""
    exp = open_expander("Choose K")
    shot(exp, "choosek.png", height=CHOOSEK_TOP)
    click("[...document.querySelectorAll('button')].find(b=>b.innerText.includes("
          "'Compare criteria against each other'))", wait=12)
    settle("document.querySelector('[role=dialog]')")
    shot("document.querySelector('[role=dialog]')", "tradeoffs.png", height=TRADEOFFS_TOP)
    click("document.querySelector('[role=dialog] button[aria-label=Close]')", wait=5)
    click(f"{exp}.querySelector('summary')", wait=4)


def report() -> None:
    """The report cards, without the shape-features card below them."""
    exp = open_expander("Report")
    shot(exp, "report.png", height=REPORT_TOP)
    click(f"{exp}.querySelector('summary')", wait=4)


def modal_shape() -> None:
    """Cell type 2's modal, shape-features view."""
    p.viewport(W, 1500, DSF); time.sleep(4)
    click(sc.LEGEND.format(2), button="right", wait=12)
    click("[...document.querySelectorAll('[role=dialog] button')].find(b=>b.innerText.includes("
          "'Shape features'))", wait=12)
    settle("document.querySelector('[role=dialog]')")
    shot("document.querySelector('[role=dialog]')", "modal_shape.png")
    click("document.querySelector('[role=dialog] button[aria-label=Close]')", wait=5)


def confidence() -> None:
    """The spatial plot coloured by mapping confidence, with the threshold slider
    and the "Plot confidence" switch under it (the UMAP beside it is unchanged)."""
    p.viewport(W, 1500, DSF); time.sleep(4)
    click("[...document.querySelectorAll('button')].find(b=>b.innerText.trim().endsWith('Plot confidence'))",
          wait=12)
    p.js("document.activeElement && document.activeElement.blur();" + NO_TIPS)
    full = sc.grab(1500)
    p.js("document.getElementById('noTip')?.remove()")
    # from the plot's title down to the switch, from the threshold label to the
    # colour bar (found on the page, so the crop follows the GUI's layout)
    title = p.rect("[...document.querySelectorAll('*')].find(e=>e.children.length==0 && "
                   "(e.textContent||'').trim()==='Spatial confidence')")
    # the threshold slider is a component drawn in a shadow DOM
    label = p.rect("(()=>{const all=[];const walk=r=>{for(const e of r.querySelectorAll('*'))"
                   "{all.push(e);if(e.shadowRoot)walk(e.shadowRoot);}};walk(document);"
                   "return all.find(e=>e.children.length==0&&(e.textContent||'').trim()"
                   ".startsWith('Confidence threshold'))})()")
    switch = p.rect("[...document.querySelectorAll('button')].find(b=>b.innerText.trim()"
                    ".endsWith('Plot confidence'))")
    x, y = label["x"] - 12, title["y"] - 12
    w, h = W - 40 - x, switch["y"] + switch["h"] + 12 - y
    full.crop((x * DSF, y * DSF, (x + w) * DSF, (y + h) * DSF)).save(OUT["confidence"] / "confidence.png")
    as_card(OUT["confidence"] / "confidence.png")
    print("wrote confidence.png")
    click("[...document.querySelectorAll('button')].find(b=>b.innerText.trim().endsWith('Plot confidence'))",
          wait=8)


def reference_markers() -> None:
    """The Single-cell reference page's marker-gene card of the Level1 annotation:
    the Level1 start at K = T, chosen on Explore (the page keeps start and K).
    Run last: it leaves the Level1 start selected."""
    sc.level1_at_t()
    click("[...document.querySelectorAll('a')].find(e=>e.innerText.trim()"
          ".endsWith('Single-cell reference'))", wait=15)
    exp = open_expander("Cell type marker genes")
    shot(exp, "ref_markers.png")


if __name__ == "__main__":
    main(set(sys.argv[1:]))
