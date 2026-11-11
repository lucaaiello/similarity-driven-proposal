SLURM MIGRATION OF THE EXISTING FRESH50K EXPERIMENT

Unchanged R/C++ sources and original input/plan hashes. New batch wrappers only.
Use a NEW remote workspace: $HOME/jcgs_dm_apt50k_slurm_v1.
Do not copy old audit receipts, native libraries, caches or production outputs.
The bundle contains original explicit source/input files plus these wrappers.
R4.4.1 on Rocky Linux; private $HOME/R/hoffman2-slurm-4.4.1 library.
Pinned Rcpp1.0.12, RcppArmadillo15.2.3.1, RcppEnsmallen0.2.22.1.1.
module purge; module load R/4.4.1. BLAS/OpenMP each limited to one thread.
Ten PSOCK processes on one node, using one Slurm task with ten allocated CPUs.
NSLOTS is exported solely for the unchanged R allocation check.

Scientific design remains50k iterations,25k burn,thin1,ten replicas,
three methods,four fresh starts,all nine adjacent edges each sweep,
original adaptive rules and finite optimization budgets. Array1-12%6;
five days=120hours and40GB per ladder. Summary1CPU/12GB/4hours.
The old Grid Engine workspace and scripts remain unchanged.

Before submission hold legacy production15351594 and summary15351595
on hoffman2.idre.ucla.edu; verify both include h in their state.
Record confirmation in slurm/LEGACY_HELD.txt in the NEW workspace.
submit.sh validates scripts and manifests, refuses existing outputs/submission
records, and queues audit->afterok ladders->afterok summary. No automatic
resubmission on failure; inspect recorded job IDs before any recovery.
The unchanged R production gate additionally requires the new runtime audit.
No launch or MCMC is performed by preparing or extracting the bundle.
Existing interactive job is not used for the batch audit and may be released.
