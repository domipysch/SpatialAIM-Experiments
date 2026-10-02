# The whole Results 2.3 figure for one mapper, from the runs on disk:
#
#   02a_analyse_invisible_split.py    the numbers    -> robustness/_analysis/<mapper>/
#   02b_panels.py                     the panel PNGs -> robustness/plots/<mapper>/, copied to
#                                                       Excalidraw/figure_3/<mapper>/panels/
#   02c_make_excalidraw.py            the drawing    -> Excalidraw/figure_3/<mapper>/figure_3_<mapper>.excalidraw.md
#
# WARNING: 02c_make_excalidraw.py overwrites the drawing, hand edits included. Run it
# for a mapper whose drawing you have since edited in Obsidian only if you mean
# to lose those edits.
#
# Every mapper gets the same figure, side by side. The runs have to exist already
# (01_run_robustness_experiments.ps1).
#
# Usage:  powershell -File 02_make_figures.ps1 -Mapper wann
#         powershell -File 02_make_figures.ps1 -All

param(
  [ValidateSet("tangram", "tacco", "dot", "nearest_centroid", "wann")]
  [string]$Mapper,
  [switch]$All
)

$ErrorActionPreference = "Continue"
if (-not $All -and -not $Mapper) { Write-Host "give -Mapper <name> or -All"; exit 2 }

$E = "C:\Users\zi69hebi\AppData\Local\miniforge3\envs\spatialaim_env"
# see 01_run_robustness_experiments.ps1: without the env on PATH the squidpy import in the
# analysis step dies with 0xC0000409
$env:PATH = "$E;$E\Library\mingw-w64\bin;$E\Library\usr\bin;$E\Library\bin;$E\Scripts;$E\bin;" +
            "C:\Users\zi69hebi\AppData\Local\miniforge3\condabin;C:\Users\zi69hebi\AppData\Local\miniforge3\Scripts;$env:PATH"

$HERE = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $HERE

$mappers = if ($All) { @("tangram", "tacco", "dot", "nearest_centroid", "wann") } else { @($Mapper) }

foreach ($m in $mappers) {
  Write-Host "`n=== $m ==="
  foreach ($step in @("02a_analyse_invisible_split.py", "02b_panels.py", "02c_make_excalidraw.py")) {
    Write-Host "--- $step"
    & "$E\python.exe" $step --mapper $m
    if ($LASTEXITCODE -ne 0) { Write-Host "  FAILED ($LASTEXITCODE), skipping $m"; break }
  }
}
