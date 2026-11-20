"""
Give a GUI screenshot the outline of a GUI card: the corners outside a rounded
rectangle on the image's bounds turn white (a modal's corners otherwise show the
dimmed page behind it), and the GUI's card border is drawn along that rectangle
(so a card cut short, or a modal, which has no border of its own, reads as one).

Radius and colour are the GUI's own (Streamlit: 0.5 rem corners, border #d6d6d9),
at the screenshots' 2 px per CSS px.
"""

from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw

RADIUS = 32                      # px: 16 CSS px
BORDER = (214, 214, 217, 255)    # the GUI's card border
WIDTH = 2                        # px: 1 CSS px


def as_card(path: Path) -> None:
    img = Image.open(path).convert("RGBA")
    box = (0, 0, img.width - 1, img.height - 1)
    inside = Image.new("L", img.size, 0)
    ImageDraw.Draw(inside).rounded_rectangle(box, RADIUS, fill=255)
    out = Image.new("RGBA", img.size, (255, 255, 255, 255))
    out.paste(img, mask=inside)
    ImageDraw.Draw(out).rounded_rectangle(box, RADIUS, outline=BORDER, width=WIDTH)
    out.convert("RGB").save(path)
