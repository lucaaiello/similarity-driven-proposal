#!/bin/bash
# Source only within a Slurm allocation. Recorded production software settings.
set -euo pipefail
: "${SLURM_JOB_ID:?Run through Slurm}"
: "${SLURM_CPUS_PER_TASK:?Explicit CPU allocation required}"
module purge
module load R/4.4.1
export R_LIBS_USER="$HOME/R/hoffman2-slurm-4.4.1"
test -d "$R_LIBS_USER"
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OMP_DYNAMIC=FALSE
export MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
# Compatibility with the allocation checks in the existing R sampler wrapper.
export NSLOTS="$SLURM_CPUS_PER_TASK"
Rscript -e 'stopifnot(getRversion()=="4.4.1"); expected<-c(Rcpp="1.0.12",RcppArmadillo="15.2.3.1",RcppEnsmallen="0.2.22.1.1"); actual<-vapply(names(expected),function(p)as.character(packageVersion(p)),""); stopifnot(identical(actual,expected));print(actual)'
sha256sum -c 'Real Data Analysis/Cluster/source_files.sha256'
