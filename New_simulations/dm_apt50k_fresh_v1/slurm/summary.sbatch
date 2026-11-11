#!/bin/bash
#SBATCH --job-name=apt50_summary
#SBATCH --account=pi_sudipto
#SBATCH --qos=highp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=04:00:00
#SBATCH --output=New_simulations/dm_apt50k_fresh_v1/slurm/logs/summary.%j.log
set -euo pipefail
source New_simulations/dm_apt50k_fresh_v1/slurm/environment.sh
Rscript New_simulations/dm_apt50k_fresh_v1/summarize.R
echo "SLURM FRESH50K SUMMARY COMPLETE job=$SLURM_JOB_ID end=$(date --iso-8601=seconds)"
