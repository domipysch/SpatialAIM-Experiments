# The three Results 2.2 sweeps on pair 1 (01_mouse-ssp x 01_mouse-ssp), all average
# linkage, for one spot-to-state mapper or all five (-All):
#
#   original_annotation          obs["cellType"] of 01_mouse-ssp, 3 types -> K = 3..1
#   invisible_split              obs["cellType_split"] of split.h5ad, 4 types -> K = 4..1
#   original_leiden_overcluster  the default Leiden over-clustering (~44 start clusters)
#
# into Data/02_Results/robustness/<run>/<mapper>/, one log per run in _logs/. A run
# whose k_selection.json exists is skipped, so an interrupted call just continues;
# delete a run's folder to redo it.
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
  @{ name = "original_annotation";         sc = $SC;       ann = "cellType" },
  @{ name = "invisible_split";             sc = $SC_SPLIT; ann = "cellType_split" },
  @{ name = "original_leiden_overcluster"; sc = $SC;       ann = $null }
)
$mappers = if ($All) { @("nearest_centroid", "wann", "tacco", "dot", "tangram") } else { @($Mapper) }

foreach ($m in $mappers) {
  foreach ($r in $runs) {
    $out = "$ROOT\$($r.name)"
    if (Test-Path "$out\$m\k_selection.json") { Write-Host "skip   $($r.name) [$m]"; continue }
    $start = if ($r.ann) { @("--start_from_annotation", $r.ann) } else { @() }
    Write-Host "run    $($r.name) [$m] ..."
    $t0 = Get-Date
    & "$E\python.exe" -m spatialaim run --scdata $r.sc --stdata $ST --output_dir $out `
        --mapping $m --linkage_method average @start *> "$ROOT\_logs\$($r.name)_$m.log"
    Write-Host ("  -> exit {0} after {1:N0} s" -f $LASTEXITCODE, ((Get-Date) - $t0).TotalSeconds)
  }
}
