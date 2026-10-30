# The three Results 2.3 sweeps on pair 1 (01_mouse-ssp x 01_mouse-ssp), all average
# linkage, for one spot-to-state mapper or all five (-All). Two spatialaim calls per
# mapper, each sweeping the Leiden start plus one annotation start:
#
#   original/annotation_cellType/<mapper>              obs["cellType"], 3 types -> K = 3..1
#   original/leiden/<mapper>                           the default Leiden over-clustering (~44 start clusters)
#   invisible_split/annotation_cellType_split/<mapper> obs["cellType_split"] of split.h5ad, 4 types -> K = 4..1
#   invisible_split/leiden/<mapper>                    by-product of the CLI (the Leiden start always runs); unused
#
# under Data/02_Results/robustness/, one log per call in _logs/. A call whose two
# k_selection.json exist is skipped, so an interrupted call just continues; delete a
# start's mapper folder to redo it.
#
# Usage:  powershell -File 01_run_robustness_experiments.ps1 -All
#         powershell -File 01_run_robustness_experiments.ps1 -Mapper wann

param(
  [ValidateSet("tangram", "tacco", "dot", "nearest_centroid", "wann")]
  [string]$Mapper,
  [switch]$All
)

$ErrorActionPreference = "Continue"
if (-not $All -and -not $Mapper) { Write-Host "give -Mapper <name> or -All"; exit 2 }

$E = "C:\Users\zi69hebi\AppData\Local\miniforge3\envs\spatialaim_env"
# Activating the env by PATH is required: calling envs\spatialaim_env\python.exe without
# it loads the wrong native DLLs and the sweep dies with 0xC0000409.
$env:PATH = "$E;$E\Library\mingw-w64\bin;$E\Library\usr\bin;$E\Library\bin;$E\Scripts;$E\bin;" +
            "C:\Users\zi69hebi\AppData\Local\miniforge3\condabin;C:\Users\zi69hebi\AppData\Local\miniforge3\Scripts;$env:PATH"

$DATA     = "C:\Users\zi69hebi\Dev\10_Alignment\Data"
$ST       = "$DATA\01_Datasets_Paper\ST\01_mouse-ssp.h5ad"
$SC       = "$DATA\01_Datasets_Paper\scRNA\01_mouse-ssp.h5ad"
$SC_SPLIT = "$DATA\01_Datasets_Paper\Robustness\split.h5ad"   # 00_perturb.py
$ROOT     = "$DATA\02_Results\robustness"
New-Item -ItemType Directory -Force -Path "$ROOT\_logs" | Out-Null
if (-not (Test-Path $SC_SPLIT)) { Write-Host "missing $SC_SPLIT - run 00_perturb.py first"; exit 1 }

$runs = @(
  @{ name = "original";        sc = $SC;       ann = "cellType" },
  @{ name = "invisible_split"; sc = $SC_SPLIT; ann = "cellType_split" }
)
$mappers = if ($All) { @("nearest_centroid", "wann", "tacco", "dot", "tangram") } else { @($Mapper) }

foreach ($m in $mappers) {
  foreach ($r in $runs) {
    $out = "$ROOT\$($r.name)"
    $done = (Test-Path "$out\leiden\$m\k_selection.json") -and
            (Test-Path "$out\annotation_$($r.ann)\$m\k_selection.json")
    if ($done) { Write-Host "skip   $($r.name) [$m]"; continue }
    Write-Host "run    $($r.name) [$m] ..."
    $t0 = Get-Date
    & "$E\python.exe" -m spatialaim run --scdata $r.sc --stdata $ST --output_dir $out `
        --mapping $m --linkage_method average --annotation $r.ann *> "$ROOT\_logs\$($r.name)_$m.log"
    Write-Host ("  -> exit {0} after {1:N0} s" -f $LASTEXITCODE, ((Get-Date) - $t0).TotalSeconds)
  }
}
