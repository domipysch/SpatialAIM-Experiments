# The whole Results 2.2 figure for one mapper, from the runs on disk:
#
#   analyse_invisible_split.py    the numbers    -> invisible_split/_analysis/<mapper>/
#   panels.py                     the panel PNGs -> Excalidraw/figure_2_2[_<mapper>]/
#   make_excalidraw.py            the drawing    -> Excalidraw/figure_2_2[_<mapper>].excalidraw.md
#   aim_experiments.render        the PDF        -> AIM-Overleaf/figures/fig_2_2[_<mapper>].pdf
#
# WARNING: make_excalidraw.py overwrites the drawing, hand edits included. Run it
# for a mapper whose drawing you have since edited in Obsidian only if you mean
# to lose those edits.
#
# Tangram keeps the unsuffixed names -- it is the figure the paper prints; every
# other mapper gets the same figure beside it. The runs have to exist already
# (run_robustness.ps1).
#
# Usage:  powershell -File make_figures.ps1 -Mapper wann
#         powershell -File make_figures.ps1 -All

param(
  [ValidateSet("tangram", "tacco", "dot", "nearest_centroid", "wann")]
  [string]$Mapper = "tangram",
  [switch]$All
)

$ErrorActionPreference = "Continue"

$E = "C:\Users\zi69hebi\AppData\Local\miniforge3\envs\aim_env"
# see run_robustness.ps1: without the env on PATH the squidpy import in the
# analysis step dies with 0xC0000409
$env:PATH = "$E;$E\Library\mingw-w64\bin;$E\Library\usr\bin;$E\Library\bin;$E\Scripts;$E\bin;" +
            "C:\Users\zi69hebi\AppData\Local\miniforge3\condabin;C:\Users\zi69hebi\AppData\Local\miniforge3\Scripts;$env:PATH"

$HERE = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $HERE

$mappers = if ($All) { @("tangram", "tacco", "dot", "nearest_centroid", "wann") } else { @($Mapper) }

foreach ($m in $mappers) {
  $suffix = if ($m -eq "tangram") { "2_2" } else { "2_2_$m" }
  Write-Host "`n=== $m ==="
  foreach ($step in @("analyse_invisible_split.py", "panels.py", "make_excalidraw.py")) {
    Write-Host "--- $step"
    & "$E\python.exe" $step --mapper $m
    if ($LASTEXITCODE -ne 0) { Write-Host "  FAILED ($LASTEXITCODE), skipping $m"; break }
  }
  if ($LASTEXITCODE -ne 0) { continue }
  Write-Host "--- render $suffix"
  & "$E\python.exe" -m aim_experiments.render $suffix
}
