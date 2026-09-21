# Both Results 2.2 experiments on pair 0 (01_mouse-ssp x 01_mouse-ssp), for one
# spot-to-state mapper. Generalises the older run_invisible_split.ps1, which did
# for Tangram alone.
#
#   invisible split   two annotation-started runs, average linkage:
#                       original  obs["cellType"]        3 types -> K = 3..1
#                       split     obs["cellType_split"]  4 types -> K = 4..1
#   granularity       the default over-clustering run (44 start clusters), plus
#                     the 3-type annotation run the older plot_granularity uses
#
# Every run is skipped when its k_selection.json is already there, so the script
# is safe to re-run and only fills what is missing.
#
# Usage:  powershell -File run_robustness.ps1 -Mapper wann
#         powershell -File run_robustness.ps1 -Mapper dot -Only granularity
#         powershell -File run_robustness.ps1 -Mapper dot -KMin 7 -KMax 7   (smoke test)

param(
  [ValidateSet("tangram", "tacco", "dot", "nearest_centroid", "wann")]
  [string]$Mapper = "tangram",
  [ValidateSet("both", "split", "granularity")]
  [string]$Only = "both",
  [int]$KMin = 0,
  [int]$KMax = 0,
  [switch]$Force
)

$ErrorActionPreference = "Continue"

$E = "C:\Users\zi69hebi\AppData\Local\miniforge3\envs\aim_env"
# Activating the env by PATH is required: calling envs\aim_env\python.exe without
# it loads the wrong native DLLs and the sweep dies with 0xC0000409.
$env:PATH = "$E;$E\Library\mingw-w64\bin;$E\Library\usr\bin;$E\Library\bin;$E\Scripts;$E\bin;" +
            "C:\Users\zi69hebi\AppData\Local\miniforge3\condabin;C:\Users\zi69hebi\AppData\Local\miniforge3\Scripts;$env:PATH"

$DATA    = "C:\Users\zi69hebi\Dev\10_Alignment\Data"
$ST      = "$DATA\01_Datasets\ST\01_mouse-ssp.h5ad"
$SC      = "$DATA\01_Datasets\scRNA\01_mouse-ssp.h5ad"
$SC_SPLIT= "$DATA\01_Datasets\Robustness\split.h5ad"
$ROOT    = "$DATA\02_Results\robustness"
$LOGS    = "$ROOT\_logs"
$LINKAGE = "average"
New-Item -ItemType Directory -Force -Path $LOGS | Out-Null

$runs = @()
if ($Only -in @("both", "split")) {
  $runs += @{ name = "invisible_split/original"; sc = $SC;       ann = "cellType" }
  $runs += @{ name = "invisible_split/split";    sc = $SC_SPLIT; ann = "cellType_split" }
}
if ($Only -in @("both", "granularity")) {
  $runs += @{ name = "granularity/annotation";  sc = $SC; ann = "cellType" }
  $runs += @{ name = "granularity/overcluster"; sc = $SC; ann = $null }
}

$ledger = "$LOGS\ledger.csv"
if (-not (Test-Path $ledger)) { "run,mapper,linkage,exit_code,seconds" | Out-File -FilePath $ledger -Encoding utf8 }

foreach ($r in $runs) {
  $outDir = "$ROOT\$($r.name)".Replace("/", "\")
  $tag    = "$($r.name)".Replace("/", "_")
  if ((-not $Force) -and (Test-Path "$outDir\$Mapper\k_selection.json")) {
    Write-Host "skip   $($r.name) [$Mapper]"; continue
  }

  $extra = @()
  if ($r.ann) { $extra += @("--start_from_annotation", $r.ann) }
  if ($KMin -gt 0) { $extra += @("--k_min", $KMin) }
  if ($KMax -gt 0) { $extra += @("--k_max", $KMax) }

  Write-Host "run    $($r.name) [$Mapper] $($extra -join ' ') ..."
  $t0 = Get-Date
  & "$E\python.exe" -m aim run `
      --scdata $r.sc --stdata $ST --output_dir $outDir `
      --mapping $Mapper --linkage_method $LINKAGE @extra `
      *> "$LOGS\${tag}_$Mapper.log"
  $code = $LASTEXITCODE
  $secs = [int]((Get-Date) - $t0).TotalSeconds
  "$($r.name),$Mapper,$LINKAGE,$code,$secs" | Out-File -FilePath $ledger -Append -Encoding utf8
  Write-Host "  -> exit $code after ${secs}s"
}

Write-Host "ROBUSTNESS RUNS DONE [$Mapper]"
foreach ($r in $runs) {
  $sel = "$ROOT\$($r.name)\$Mapper\k_selection.json".Replace("/", "\")
  if (Test-Path $sel) {
    $best = (Get-Content $sel -Raw | ConvertFrom-Json).best.overall.k
    Write-Host ("  {0,-28} best overall K = {1}" -f $r.name, $best)
  }
}
