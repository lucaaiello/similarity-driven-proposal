#!/bin/bash
#SBATCH --job-name=dm_analysis
#SBATCH --account=pi_sudipto
#SBATCH --qos=highp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --time=5-00:00:00
#SBATCH --array=1-12%6
#SBATCH --output="Real Data Analysis/Cluster/Logs/analysis.%A.%a.log"
set -euo pipefail
source 'Real Data Analysis/Cluster/environment.sh'
test "$SLURM_CPUS_PER_TASK" -eq 10
echo "ANALYSIS job=$SLURM_JOB_ID task=$SLURM_ARRAY_TASK_ID host=$(hostname) start=$(date --iso-8601=seconds)"
Rscript 'Real Data Analysis/run_analysis.R' --task="$SLURM_ARRAY_TASK_ID" --workers="$SLURM_CPUS_PER_TASK"
echo "ANALYSIS COMPLETE job=$SLURM_JOB_ID task=$SLURM_ARRAY_TASK_ID end=$(date --iso-8601=seconds)"
