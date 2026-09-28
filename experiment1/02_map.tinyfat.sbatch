#!/bin/bash -l
# Stage B, one mapper (nearest_centroid | wann | tacco | dot | tangram, all on CPU) for one unit:
# ward, then average. Submit from this directory, once stage A is complete:
#   sbatch.tinyfat -J exp1_B_wann --time=02:00:00 --array=0-251%20 02_map.tinyfat.sbatch run wann
# --time per mapper and start kind: grid.job_time (the default below is 2 h, DOT gets more).
#SBATCH --job-name=exp1_B_cpu
#SBATCH --partition=work
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
bash -l "$SLURM_SUBMIT_DIR/run_job.sh" "${1:-run}" B "${2:?mapper}"
