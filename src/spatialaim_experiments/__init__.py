"""Helpers shared by the figure and experiment scripts of the SpatialAIM paper.

    paths        every filesystem location, with environment overrides
    style        palette, mapper names, and PanelWriter (panel PNG + manifest)
    excalidraw   building the drawing each figure is assembled as
    exp1         the Experiment 1 results table (pulled summary)

The SpatialAIM method itself is not here: it is the ``spatial-aim`` package (the
SpatialAIM-Repository), which this repository imports as ``spatialaim`` and never
modifies.
"""

__all__ = ["excalidraw", "exp1", "paths", "style"]
