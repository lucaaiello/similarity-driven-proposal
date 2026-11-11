#!/bin/bash
#SBATCH --job-name=apt50_audit
#SBATCH --account=pi_sudipto
#SBATCH --qos=highp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --time=03:00:00
#SBATCH --output=New_simulations/dm_apt50k_fresh_v1/slurm/logs/audit.%j.log
set -euo pipefail
source New_simulations/dm_apt50k_fresh_v1/slurm/environment.sh
test "$SLURM_CPUS_PER_TASK" -eq 10
echo "SLURM AUDIT job=$SLURM_JOB_ID host=$(hostname) cores=$SLURM_CPUS_PER_TASK start=$(date --iso-8601=seconds)"
Rscript New_simulations/dm_apt_runtime_optimization/verify.R --mode=full --workers=10 --output=New_simulations/dm_apt50k_fresh_v1/engine_audit
Rscript New_simulations/dm_apt50k_fresh_v1/audit.R --workers=10 --output=New_simulations/dm_apt50k_fresh_v1/audit
echo "SLURM FRESH50K AUDITS COMPLETE job=$SLURM_JOB_ID end=$(date --iso-8601=seconds)"
