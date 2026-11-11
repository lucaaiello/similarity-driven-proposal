# Diet-microbiome analysis

This folder contains the application data, final sampler, saved results,
diagnostics, figure/table code and optional cluster workflow.

```text
Real Data Analysis/
  Data/                     # supplied counts/covariates/labels
    Prepared/               # exact fresh inputs and analysis plan
    Reference/              # states used only in implementation checks
  Sampler/                  # native sampler and parallel-tempering wrapper
  Results/                  # completed 12-ladder cold-chain summary and CSVs
    Diagnostics/            # saved extended assessment
    Table_Summaries/        # saved manuscript numerical exports
  Figures_and_Tables/       # saved-results plotting and table construction
  Validation/               # checks and independent/reference implementations
  Cluster/                  # Slurm wrappers and source/input bundle preparation
  Provenance/               # original source archive, receipts and location map
  run_analysis.R            # one new sampling task, only when explicitly invoked
  summarize_results.R       # summarize new complete checkpoints
  assess_results.R          # extra local diagnostics from the saved summary
  prepare_inputs.R          # pure-R verification of the supplied fresh starts
```

## Reproduce the reported results

From the repository root, use:

```sh
Rscript reproduce/run.R
```

This also reproduces the simulation displays and verifies every computational
table. It runs no MCMC, optimizer, compiler or cluster job. For application-only
assets in a new destination:

```sh
Rscript "Real Data Analysis/Figures_and_Tables/create_figures.R" --output-root=reproduced_application
```

`Results/summary.rds` contains the completed results for three proposal methods,
four independent ladders per method, ten replicas per ladder, 50,000 sweeps per
replica, burn-in 25,000, and no thinning. Each ladder contributes one retained
cold posterior chain. All 30 genera are included in the fitted model. The
without-Prevotella summaries describe the same full-model draws, not refits.

The data and numerical results are unchanged. Personal paths in saved diagnostic
metadata have been redacted; no numerical values were changed. Original source
receipts are checked against `Provenance/original_sources.zip`, whose archived
files retain their recorded names. Current live source paths are deliberately
separate and require fresh implementation checks before a new sampler run.
This separation preserves historical provenance without pretending that a
renamed executable file has its old checksum. No result receipt was rewritten.

## Run the application samplers again

These commands generate new MCMC results and can take days on the cluster.
They are not needed to rebuild tables and figures from the supplied results.
Read `Cluster/README.md` and `Validation/README.md` first. New runs write to
`Results/New_Run/`; they do not overwrite `Results/summary.rds`. The supplied
prepared inputs define the original fresh null/5%/10%/15% starts. Running
`prepare_inputs.R` verifies their reconstruction without changing the files.

Full remote continuation checkpoints were not all downloaded. The saved summary
supports every reported table/figure but is not a complete restart state for
each hot replica. The original serialized schema identifiers remain inside
records for compatibility; active filenames do not encode iteration counts.
