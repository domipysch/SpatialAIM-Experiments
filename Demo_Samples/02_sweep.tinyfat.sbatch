#!/bin/bash -l
# One run of runs.tsv (pair x start x linkage x mapper) per array task: the full K sweep.
# Submitted per walltime class by 03_submit_sweeps.sh, e.g.
#   sbatch.tinyfat -J demo_tangram --time=02:00:00 --array=$(python grid.py array tangram)%20 02_sweep.tinyfat.sbatch
#SBATCH --job-name=demo_sweep
#SBATCH --partition=work
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
bash -l "$SLURM_SUBMIT_DIR/run_job.sh" sweep "$SLURM_ARRAY_TASK_ID"
