#!/bin/bash -l
# Submit every sweep, once all scaffolds are done (python run_one.py status): one array per
# walltime class, plus one per-K array for each run split by K. Run from this directory with
# the env active (source config.sh && conda activate $SPATIALAIM_ENV).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
mkdir -p logs

for cls in $(python -c "import grid; print(*grid.CLASSES)"); do
  array="$(python grid.py array "$cls")"
  [[ -n "$array" ]] || continue
  sbatch.tinyfat -J "demo_$cls" --time="$(python grid.py time "$cls")" --array="${array}%20" 02_sweep.tinyfat.sbatch
done

for run in $(python -c "import grid; print(*[r['run'] for r in grid.runs() if grid.split_by_k(r)])"); do
  sbatch.tinyfat --array="$(python run_one.py levels "$run")%20" 02b_sweep_by_k.tinyfat.sbatch "$run"
done
