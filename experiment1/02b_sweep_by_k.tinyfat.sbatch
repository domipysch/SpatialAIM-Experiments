#!/bin/bash -l
# One K of one sweep per array task (the array index is K), for a sweep that does not
# fit a 24 h job (so far: DOT on pair 93, 09_human-breast-cancer, Leiden start).
# Submit from this directory, per tree without DONE:
#   sbatch.tinyfat --array=$(python run_one.py levels <unit>)%20 02b_sweep_by_k.tinyfat.sbatch run <unit> <mapper> <tree>
# then, once every K is done: python run_one.py merge <unit> <mapper> <tree>
#SBATCH --job-name=exp1_sweep_k
#SBATCH --partition=work
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
bash -l "$SLURM_SUBMIT_DIR/run_job.sh" "${1:-run}" k "${2:?unit}" "${3:?mapper}" "${4:?tree}" "$SLURM_ARRAY_TASK_ID"
