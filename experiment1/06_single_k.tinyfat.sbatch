#!/bin/bash -l
# Runtime of mapping one annotation (K = T only) for one annotation unit and mapper;
# the array index is the unit. Submit from this directory, once per mapper:
#   sbatch.tinyfat -J exp1_T_wann --array=$(python grid.py array annotation)%20 06_single_k.tinyfat.sbatch run wann
# Same node type, cores, threads and memory as stage B (02_map.tinyfat.sbatch), so the times compare.
#SBATCH --job-name=exp1_T
#SBATCH --partition=work
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
source "$SLURM_SUBMIT_DIR/config.sh" || exit 1
case "${1:?root}" in
  run|test) export EXP1_ROOT="$EXP1_BASE/$1" ;;
  *)        export EXP1_ROOT="$1" ;;
esac
conda activate "$SPATIALAIM_ENV"
echo "=== $(date -Is)  job $SLURM_JOB_ID (array ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID})  $(hostname)"
python "$SLURM_SUBMIT_DIR/06_single_k.py" "$SLURM_ARRAY_TASK_ID" "${2:?mapper}"
