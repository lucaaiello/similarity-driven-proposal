#!/bin/bash
#SBATCH --job-name=apt50_fresh
#SBATCH --account=pi_sudipto
#SBATCH --qos=highp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --time=5-00:00:00
#SBATCH --array=1-12%6
#SBATCH --output=New_simulations/dm_apt50k_fresh_v1/slurm/logs/ladder.%A.%a.log
set -euo pipefail
source New_simulations/dm_apt50k_fresh_v1/slurm/environment.sh
test "$SLURM_CPUS_PER_TASK" -eq 10
echo "SLURM FRESH50K job=$SLURM_JOB_ID task=$SLURM_ARRAY_TASK_ID host=$(hostname) cores=$SLURM_CPUS_PER_TASK start=$(date --iso-8601=seconds)"
Rscript New_simulations/dm_apt50k_fresh_v1/run_task.R --task="$SLURM_ARRAY_TASK_ID" --workers="$SLURM_CPUS_PER_TASK"
echo "SLURM FRESH50K COMPLETE job=$SLURM_JOB_ID task=$SLURM_ARRAY_TASK_ID end=$(date --iso-8601=seconds)"
