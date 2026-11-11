#!/bin/bash
#SBATCH --job-name=dm_audit
#SBATCH --account=pi_sudipto
#SBATCH --qos=highp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --time=03:00:00
#SBATCH --output="Real Data Analysis/Cluster/Logs/audit.%j.log"
set -euo pipefail
source 'Real Data Analysis/Cluster/environment.sh'
test "$SLURM_CPUS_PER_TASK" -eq 10
echo "AUDIT job=$SLURM_JOB_ID host=$(hostname) start=$(date --iso-8601=seconds)"
Rscript 'Real Data Analysis/Validation/check_engine.R' --mode=full --workers=10 --output='Real Data Analysis/Validation/Results/Engine'
Rscript 'Real Data Analysis/Validation/check_analysis.R' --workers=10 --output='Real Data Analysis/Validation/Results/Analysis'
echo "AUDITS COMPLETE job=$SLURM_JOB_ID end=$(date --iso-8601=seconds)"
