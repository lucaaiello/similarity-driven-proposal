#!/bin/bash
#SBATCH --job-name=dm_summary
#SBATCH --account=pi_sudipto
#SBATCH --qos=highp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=04:00:00
#SBATCH --output="Real Data Analysis/Cluster/Logs/summary.%j.log"
set -euo pipefail
source 'Real Data Analysis/Cluster/environment.sh'
Rscript 'Real Data Analysis/summarize_results.R'
echo "SUMMARY COMPLETE job=$SLURM_JOB_ID end=$(date --iso-8601=seconds)"
