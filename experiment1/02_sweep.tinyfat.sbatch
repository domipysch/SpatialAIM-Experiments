#!/bin/bash -l
# One mapper (nearest_centroid | wann | tacco | dot | tangram, all on CPU) for one unit:
# the full K sweep, ward then average. Submit from this directory, once the
# scaffolds are complete:
#   sbatch.tinyfat -J exp1_sweep_wann --time=02:00:00 --array=0-192%20 02_sweep.tinyfat.sbatch run wann
# --time per mapper and start kind: grid.job_time (the default below is 2 h, DOT gets more).
#SBATCH --job-name=exp1_sweep
#SBATCH --partition=work
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
bash -l "$SLURM_SUBMIT_DIR/run_job.sh" "${1:-run}" sweep "$SLURM_ARRAY_TASK_ID" "${2:?mapper}"
