"""
GUI screenshots of figure 5, taken at 2x device scale from a headless Chrome.

Needs the GUI running (``spatialaim gui --server_port 8503``) and Chrome started
with remote debugging:

    chrome --headless=new --remote-debugging-port=9333 --user-data-dir=<tmp> about:blank

Writes to Excalidraw/figure_5/screenshots/ (CSS-pixel crops, 2 px per CSS px):

    explore_k10.png   (a) Explore page, Leiden start, K = 10
    modal_spp1.png    (c) cell-type modal of cell type 2, SPP1; modal_spp1_umap.png its
                      reference-UMAP map
    modal_stab1.png   (c) the same for cell type 3, STAB1; modal_stab1_umap.png
    slider_k16.png    (f) the K controls at K = 16, "next Pareto K" tooltip showing
    umap_k10.png      (f) reference UMAP, K = 10
    umap_k16.png      (f) reference UMAP, K = 16
    umap_level1.png   (f) reference UMAP, Level1 annotation (K = 16 colouring)
    spatial_k10.png, spatial_k16.png, spatial_level1.png
                      (f) the same three on the Xenium section (mapped cell types)

and callouts.json: the CSS positions of what (a)'s numbered callouts point at.
Supplementary figure S7's screenshots: screenshots_s7.py, which reuses these helpers.
"""

from __future__ import annotations

import io
import json
import sys
import time
from pathlib import Path

from PIL import Image

sys.path.insert(0, str(Path(__file__).parent))
from cdp import Page  # noqa: E402
from spatialaim_experiments import paths  # noqa: E402

URL = "http://localhost:8503"
GUI = paths.DATA / "03_GUI"
OUT = paths.figure_dir("5") / "screenshots"
DSF = 2
W = 2000
FIELDS = {
    "scRNA-seq .h5ad path": GUI / "scRNA" / "11_human-crc-p2.h5ad",
    "ST .h5ad path": GUI / "ST" / "11_human-crc-p2.h5ad",
    "Segmentation path (optional)": GUI / "ST" / "Segmentation" / "11_human-crc-p2.tif",
    "Output directory": GUI / "Result",
}
# crops in CSS px (x, y, w, h)
SLIDER_CROP = (540, 95, 1400, 150)            # slider + Pareto buttons (old layout)
K_SLIDER_CROP = (540, 125, 1040, 105)         # the K slider alone, label to end ticks
UMAP_CROP = (150, 430, 690, 630)
SPATIAL_CROP = (1160, 430, 620, 650)          # "Mapped cell types" on the Xenium section
# the annotation row's two plots get the same crops, shifted down to that row
COMPUTED_UMAP_TITLE = "Reference UMAP — computed"
COMPUTED_SPATIAL_TITLE = "Mapped cell types"
LEVEL1_UMAP_TITLE = "Reference UMAP — Level1"
LEVEL1_SPATIAL_TITLE = "Mapped Level1"
FULL = paths.outputs("05_gui")                # the uncropped pages, for re-cropping
MODAL_UMAP_CROP = (288, 493, 485, 460)   # relative to the dialog: the gene on the reference UMAP
MODAL_GENES = [(2, "SPP1"), (3, "STAB1")]   # K = 10: SPP1+ and SELENOP+ macrophages
# (a)'s numbered callouts, in reading order: the element each sits beside (JS),
# found on the page so they follow the GUI's layout
CALLOUT_TARGETS = [
    "[...document.querySelectorAll('[data-testid=stSidebar] h3, [data-testid=stSidebar] p')].find(e=>e.innerText.trim()==='1. Select data')",
    "[...document.querySelectorAll('[data-testid=stSidebar] h3, [data-testid=stSidebar] p')].find(e=>e.innerText.trim()==='2. Configure method')",
    "[...document.querySelectorAll('[data-testid=stMain] p')].find(e=>e.innerText.trim()==='Starting clusters')",
    "(()=>{const all=[];const walk=r=>{for(const e of r.querySelectorAll('*')){all.push(e);if(e.shadowRoot)walk(e.shadowRoot);}};walk(document);for(const f of document.querySelectorAll('iframe')){try{walk(f.contentDocument)}catch(_){}};return all.find(e=>e.children.length==0&&(e.textContent||'').trim().startsWith('K (number'))})()",
    "[...document.querySelectorAll('button')].find(b=>b.innerText.trim()==='← Pareto K')",
    "[...document.querySelectorAll('[data-testid=stExpander] summary')].find(s=>s.innerText.includes('Choose K'))",
    "document.querySelector('.legendtitletext')",
    "[...document.querySelectorAll('*')].find(e=>e.children.length==0 && (e.textContent||'').startsWith('Reference UMAP — computed'))",
    "[...document.querySelectorAll('*')].find(e=>e.children.length==0 && (e.textContent||'').startsWith('Mapped cell types'))",
]
# (a) is cropped from below the top navigation bar to just above the plots' hint
# line and toolbar
TOP_CUT = 80
BOTTOM_CUT = 115  # CSS px above the toolbar: drops the hint line, axis titles and x ticks
TOOLBAR = "[...document.querySelectorAll('button')].find(b=>b.innerText.trim().endsWith('Shared-gene UMAP'))"

p = Page()


def click(js: str, button: str = "left", wait: float = 10) -> None:
    r = p.rect(js)
    if r is None:
        raise RuntimeError(f"not found: {js}")
    p.click(r["x"] + r["w"] / 2, r["y"] + r["h"] / 2, button=button)
    time.sleep(wait)


def grab(height: int) -> Image.Image:
    p.viewport(W, height, DSF)
    time.sleep(5)
    p.move(W - 5, height - 5)
    time.sleep(2)
    tmp = OUT / "_full.png"
    p.shot(str(tmp), clip=(0, 0, W, height))
    img = Image.open(io.BytesIO(tmp.read_bytes()))
    tmp.unlink()
    return img


def crop(img: Image.Image, box, name: str) -> None:
    x, y, w, h = box
    img.crop((x * DSF, y * DSF, (x + w) * DSF, (y + h) * DSF)).save(OUT / name)
    print("wrote", name)


def status() -> str:
    """The status line under the controls, e.g. "Leiden · nearest_centroid · K = 10 ..."."""
    return p.js("[...document.querySelectorAll('p,div')].map(e=>e.innerText)"
                ".find(t=>t && /^[\\w ]+ · nearest_centroid · K = /.test(t))") or ""


def title_y(text: str) -> float:
    r = p.rect(f"[...document.querySelectorAll('*')].find(e=>e.children.length==0 && "
               f"(e.textContent||'').startsWith('{text}'))")
    if r is None:
        raise RuntimeError(f"title not found: {text}")
    return r["y"]


def shifted(box, dy: float):
    x, y, w, h = box
    return (x, int(y + dy), w, h)


LEGEND = "[...document.querySelectorAll('.legendtext')].find(e=>e.textContent.trim()==='Cell type {}')"


def toggle_states(states) -> None:
    """Left click legend entries: shows / hides those cell types in both plots."""
    for s in states:
        click(LEGEND.format(s), wait=2)
    time.sleep(4)


def modal_gene(state: int, gene: str) -> None:
    """Open cell type ``state``'s modal, select ``gene`` from whichever of its
    marker lists carries it, and save the modal plus its reference-UMAP map."""
    p.viewport(W, 1500, DSF); time.sleep(4)
    click(f"[...document.querySelectorAll('.legendtext')].find(e=>e.textContent.trim()==='Cell type {state}')",
          button="right", wait=12)
    button = f"[...document.querySelectorAll('[role=dialog] button')].find(b=>b.innerText.trim().startsWith('{gene} '))"
    for listing in ("Upregulated · ST", "Upregulated · scRNA-seq"):
        click("document.querySelector('[role=dialog] [data-testid=stSelectbox] [role=combobox]')", wait=2)
        click(f"[...document.querySelectorAll('[role=option]')].find(o=>o.innerText.includes('{listing}'))", wait=12)
        if p.rect(button):
            break
    else:
        raise RuntimeError(f"{gene} is in none of cell type {state}'s marker lists")
    click(button, wait=12)
    # drop the button's hover/focus tooltip before the shot
    d = p.rect("document.querySelector('[role=dialog]')")
    p.js("document.activeElement && document.activeElement.blur();"
         "const s=document.createElement('style'); s.id='noTip';"
         "s.textContent='[role=tooltip],[data-baseweb=tooltip],[data-baseweb=popover]{display:none!important}';"
         "document.head.appendChild(s)")
    p.move(d["x"] + 300, d["y"] + 30); time.sleep(3)
    name = f"modal_{gene.lower()}"
    p.shot(str(OUT / f"{name}.png"), clip=(d["x"], d["y"], d["w"], d["h"]))
    print(f"wrote {name}.png ({listing})")
    crop(Image.open(OUT / f"{name}.png"), MODAL_UMAP_CROP, f"{name}_umap.png")
    p.js("document.getElementById('noTip')?.remove()")
    click("document.querySelector('[role=dialog] button[aria-label=Close]')", wait=5)


def expander_shot(title: str, name: str, height: int = 3600, out: Path = OUT) -> None:
    """Open the Explore page's expander whose header starts with ``title``, wait
    for its lazily computed content, and save the open expander."""
    p.viewport(W, height, DSF); time.sleep(4)
    exp = (f"[...document.querySelectorAll('[data-testid=stExpander]')]"
           f".find(e=>{{const s=e.querySelector('summary'); return s && s.innerText.includes('{title}')}})")
    click(f"{exp}.querySelector('summary')", wait=8)
    for _ in range(60):   # marker rankings run on a background thread
        if not p.js(f"{exp}.innerText.includes('Computing')"):
            break
        time.sleep(3)
    time.sleep(8)
    p.move(W - 5, 5); time.sleep(2)
    r = p.rect(exp)
    p.shot(str(out / name), clip=(r["x"], r["y"], r["w"], min(r["h"], height - r["y"])))
    print(f"wrote {name} ({r['w']:.0f} x {r['h']:.0f} CSS px)")
    click(f"{exp}.querySelector('summary')", wait=4)


def open_result() -> None:
    """Load the CRC P2 run through the sidebar and wait for the Explore page at
    the Leiden start's best K = 10."""
    OUT.mkdir(parents=True, exist_ok=True)
    p.viewport(1600, 1000, DSF)
    p.goto(URL, wait=8)
    # The sidebar remembers being collapsed; open it to reach the path fields.
    r = p.rect('document.querySelector(\'input[aria-label="Output directory"]\')')
    if r is None or r["w"] == 0 or r["x"] < 0:
        click("document.querySelector('[data-testid=\"stExpandSidebarButton\"]') || "
              "document.querySelector('[data-testid=\"stSidebarCollapsedControl\"] button')", wait=3)
    for label, value in FIELDS.items():
        r = p.rect(f'document.querySelector(\'input[aria-label="{label}"]\')')
        p.click(r["x"] + r["w"] / 2, r["y"] + r["h"] / 2, count=3)
        p.type(str(value)); p.key("Enter"); time.sleep(3)
    for _ in range(60):  # loading the results takes a while
        if "K = 10" in status():
            break
        time.sleep(3)
    time.sleep(8)
    assert "K = 10" in status(), status()


def collapse_sidebar() -> None:
    """Its button only reacts while hovered."""
    expanded = "document.querySelector('[data-testid=stSidebar]').getAttribute('aria-expanded')"
    for _ in range(3):
        if p.js(expanded) == "false":
            break
        p.move(150, 300); time.sleep(1)
        click("document.querySelector('[data-testid=\"stSidebarCollapseButton\"] button')", wait=4)
    assert p.js(expanded) == "false", "sidebar did not collapse"


def level1_at_t() -> None:
    """On the Explore page: switch to the Level1 start and step to its K = T = 9,
    where the tree cut is the identity (the cell types are the Level1 labels)."""
    p.move(5, 5)
    click("document.querySelector('[data-testid=stMain] [data-testid=stSelectbox] [role=combobox]')", wait=2)
    click("[...document.querySelectorAll('[role=option]')].find(o=>o.innerText.trim().endsWith('Level1'))", wait=12)
    for _ in range(10):
        if "K = 9" in status():
            break
        r = p.rect("[...document.querySelectorAll('button')].find(b=>b.innerText.trim()==='Pareto K →')")
        if r is None or p.js("[...document.querySelectorAll('button')].find(b=>b.innerText.trim()==='Pareto K →').disabled"):
            break
        p.click(r["x"] + r["w"] / 2, r["y"] + r["h"] / 2); time.sleep(10)
    print("Level1 start:", status())
    assert "K = 9" in status(), status()


def main() -> None:
    open_result()

    # (a) the Explore page, sidebar open, down to the plots
    full = grab(1500)
    toolbar = p.rect(TOOLBAR)
    explore_crop = (0, TOP_CUT, W, int(toolbar["y"] - BOTTOM_CUT - TOP_CUT))
    crop(full, explore_crop, "explore_k10.png")
    points = {}
    for n, js in enumerate(CALLOUT_TARGETS, start=1):
        r = p.rect(js)
        points[str(n)] = (max(r["x"] - 26, 16), r["y"] + r["h"] / 2)
    (OUT / "callouts.json").write_text(json.dumps({"crop": explore_crop, "points": points}))
    collapse_sidebar()   # for the rest

    # (c) the Explore UMAP with only the two macrophage states shown
    p.viewport(W, 1500, DSF); time.sleep(4)
    others = [s for s in range(10) if s not in (2, 3)]
    toggle_states(others)
    crop(grab(1500), UMAP_CROP, "umap_k10_states23.png")
    p.viewport(W, 1500, DSF); time.sleep(3)
    toggle_states(others)

    # (c) the SPP1+ marker in cell type 2's modal, the SELENOP+ marker in cell type 3's
    for state, gene in MODAL_GENES:
        modal_gene(state, gene)

    # (e) annotation row on; K = 10, then K = 16 via the Pareto button
    click("[...document.querySelectorAll('button')].find(b=>b.innerText.trim().endsWith('Show annotation'))", wait=12)
    pages = {}
    for k in (10, 16):
        p.viewport(W, 1500, DSF); time.sleep(4)
        for _ in range(6):
            if f"K = {k}," in status() or f"K = {k} " in status():
                break
            click("[...document.querySelectorAll('button')].find(b=>b.innerText.trim()==='Pareto K →')")
        full = grab(2500)
        full.save(FULL / f"page_k{k}_annotation.png")
        crop(full, UMAP_CROP, f"umap_k{k}.png")
        crop(full, SPATIAL_CROP, f"spatial_k{k}.png")
        crop(full, K_SLIDER_CROP, f"kslider_k{k}.png")
        pages[k] = full
    # the static Level1 annotation: its row of the K = 16 page (colours matched to
    # K = 16), cut exactly like the computed plots above it
    dy_u = title_y(LEVEL1_UMAP_TITLE) - title_y(COMPUTED_UMAP_TITLE)
    dy_s = title_y(LEVEL1_SPATIAL_TITLE) - title_y(COMPUTED_SPATIAL_TITLE)
    crop(pages[16], shifted(UMAP_CROP, dy_u), "umap_level1.png")
    crop(pages[16], shifted(SPATIAL_CROP, dy_s), "spatial_level1.png")

    # the tooltip of the next-Pareto button, as a hint that the controls are live
    p.viewport(W, 1500, DSF); time.sleep(4)
    r = p.rect("[...document.querySelectorAll('button')].find(b=>b.innerText.trim()==='Pareto K →')")
    p.move(r["x"] + r["w"] / 2, r["y"] + r["h"] / 2); time.sleep(3)
    tmp = OUT / "_full.png"
    p.shot(str(tmp), clip=(0, 0, W, 1500))
    crop(Image.open(io.BytesIO(tmp.read_bytes())), SLIDER_CROP, "slider_k16.png")
    tmp.unlink()

    # the Level1 start's slider at K = T = 9, whose mapping is the annotation itself
    level1_at_t()
    p.move(5, 5); time.sleep(2)
    tmp = OUT / "_full.png"
    p.shot(str(tmp), clip=(0, 0, W, 1500))
    crop(Image.open(io.BytesIO(tmp.read_bytes())), K_SLIDER_CROP, "kslider_level1.png")
    tmp.unlink()


if __name__ == "__main__":
    main()
