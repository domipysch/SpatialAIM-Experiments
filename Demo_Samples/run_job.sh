#!/bin/bash -l
# Body of every demo-results array task: set up the env, then run_one.py with the arguments.
#
#   run_job.sh scaffold <pair_id>        (called by the *.sbatch scripts)
#   run_job.sh sweep <run>
#   run_job.sh k <run> <K>
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "$HERE/config.sh" || exit 1
conda activate "$SPATIALAIM_ENV"
set -uo pipefail
mkdir -p "$DEMO_OUT" "$DEMO_STATE"

echo "=== $(date -Is)  job $SLURM_JOB_ID (array ${SLURM_ARRAY_JOB_ID:-}_${SLURM_ARRAY_TASK_ID:-})  $(hostname)"
echo "$*   cpus ${SLURM_CPUS_PER_TASK:-?}   threads $DEMO_THREADS   TMPDIR ${TMPDIR:-}"
python "$HERE/run_one.py" "$@"
