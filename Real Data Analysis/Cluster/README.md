# Optional Slurm production workflow

Nothing in this folder is invoked by saved-results reproduction. The shell
wrappers submit jobs only when `submit.sh` is explicitly run on the cluster.

The design and resource requests retain the completed application's settings:
12 ladders (three methods x four fresh starts), ten one-thread workers per
ladder, concurrency limit six, 40 GB and 120 hours per ladder. The job array
interleaves the three methods within each start. Account `pi_sudipto`, highp QOS,
R 4.4.1 and the private library path are recorded Hoffman2 settings, not generic
settings guaranteed for another cluster. Adjust wrappers for your account/site.

## Prepare files locally

From the repository root:

```sh
python "Real Data Analysis/Cluster/prepare_bundle.py"
```

This writes a verified source/input archive and checksum under `dist/`; no R,
compiler, MCMC, network access or submission is performed. Source/input files
are hashed in `source_files.sha256`. Large saved results and checkpoints are
excluded. Every archive file has owner-write permissions for safe extraction.

## New compute workspace

Extract the bundle into a NEW cluster directory. Load R 4.4.1 and set
`R_LIBS_USER="$HOME/R/hoffman2-slurm-4.4.1"`. If needed, install the pinned native
packages with `install_packages.R` inside a compute allocation. The environment
wrapper limits worker BLAS/OpenMP threads to one. Rcpp/RcppArmadillo/
RcppEnsmallen versions are checked before execution.

From that workspace root, explicitly submit with:

```sh
bash "Real Data Analysis/Cluster/submit.sh"
```

The chain is full engine/analysis audits -> afterok array -> afterok summary.
Every task also enforces matching-runtime/source audit receipts. Existing
submission records, check outputs or run checkpoints stop duplicate submission.
Logs and submitted job IDs remain under `Cluster/`; new production states and
summary go under `Results/New_Run/`. The original completed saved summary is not
overwritten. Do not resubmit merely because jobs are pending.

Local reorganization verification covers source parsing, unchanged native
bytes, original input/source receipts and saved-results reproduction. It does
not execute the full MCMC audit or claim a newly tested Slurm deployment. Those
checks run in the target cluster runtime before any new production task.
