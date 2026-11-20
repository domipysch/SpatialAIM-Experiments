"""Minimal Chrome DevTools Protocol driver for taking high-resolution GUI
screenshots of figure 5 (Chrome must run with --remote-debugging-port=9333,
headless, --force-device-scale-factor=2).

    from cdp import Page
    p = Page()                       # first ordinary tab
    p.goto("http://localhost:8503")
    p.click(x, y); p.type("text"); p.key("Enter")
    p.shot("out.png", clip=(x, y, w, h))   # CSS px; saved at device scale
"""

from __future__ import annotations

import base64
import itertools
import json
import time
import urllib.request

from websockets.sync.client import connect

PORT = 9333


class Page:
    def __init__(self, port: int = PORT):
        tabs = json.load(urllib.request.urlopen(f"http://127.0.0.1:{port}/json"))
        tab = next(t for t in tabs if t["type"] == "page")
        self.ws = connect(tab["webSocketDebuggerUrl"], max_size=None)
        self.ids = itertools.count(1)

    def send(self, method: str, **params):
        i = next(self.ids)
        self.ws.send(json.dumps({"id": i, "method": method, "params": params}))
        while True:
            msg = json.loads(self.ws.recv())
            if msg.get("id") == i:
                if "error" in msg:
                    raise RuntimeError(msg["error"])
                return msg.get("result", {})

    def viewport(self, w: int, h: int, scale: float = 2) -> None:
        self.send("Emulation.setDeviceMetricsOverride", width=w, height=h,
                  deviceScaleFactor=scale, mobile=False)

    def goto(self, url: str, wait: float = 5) -> None:
        self.send("Page.navigate", url=url)
        time.sleep(wait)

    def js(self, expr: str):
        r = self.send("Runtime.evaluate", expression=expr, returnByValue=True, awaitPromise=True)
        return r.get("result", {}).get("value")

    def click(self, x: float, y: float, button: str = "left", count: int = 1) -> None:
        for kind in ("mouseMoved", "mousePressed", "mouseReleased"):
            self.send("Input.dispatchMouseEvent", type=kind, x=x, y=y,
                      button="none" if kind == "mouseMoved" else button,
                      buttons=0 if kind == "mouseMoved" else (1 if button == "left" else 2),
                      clickCount=0 if kind == "mouseMoved" else count)

    def move(self, x: float, y: float) -> None:
        self.send("Input.dispatchMouseEvent", type="mouseMoved", x=x, y=y)

    def type(self, text: str) -> None:
        self.send("Input.insertText", text=text)

    def key(self, name: str) -> None:
        codes = {"Enter": 13, "Escape": 27, "Tab": 9}
        for kind in ("keyDown", "keyUp"):
            self.send("Input.dispatchKeyEvent", type=kind, key=name, code=name,
                      windowsVirtualKeyCode=codes[name], text="\r" if name == "Enter" and kind == "keyDown" else "")

    def shot(self, path: str, clip: tuple[float, float, float, float] | None = None,
             scale: float = 1) -> None:
        params = {"format": "png", "captureBeyondViewport": False}
        if clip is not None:
            x, y, w, h = clip
            params["clip"] = {"x": x, "y": y, "width": w, "height": h, "scale": scale}
        data = self.send("Page.captureScreenshot", **params)["data"]
        with open(path, "wb") as fh:
            fh.write(base64.b64decode(data))

    def rect(self, selector_js: str):
        """Bounding box {x,y,w,h} of the element a JS expression returns."""
        return self.js(f"(()=>{{const e={selector_js}; if(!e) return null;"
                       "const r=e.getBoundingClientRect(); return {x:r.x,y:r.y,w:r.width,h:r.height};})()")
