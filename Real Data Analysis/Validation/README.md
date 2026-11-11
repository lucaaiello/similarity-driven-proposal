# Application implementation checks

`check_engine.R` compares the native engine with an independent/reference
implementation. `check_analysis.R` checks the fresh-start, exchange, parallel RNG
and serialized-restart machinery. Reference implementations and deterministic
fixtures are contained here rather than mixed with the production entry point.

The checks have distinct scopes:

- `Rscript "Real Data Analysis/Validation/check_engine.R" --mode=parse` checks
  executable source parsing only, without compilation or MCMC.
- `--mode=compile` compiles the engines; `fixed` uses bounded fitting/helper
  comparisons; `quick`/`full` additionally run bounded synthetic chains.
- `check_analysis.R` performs finite-state checks and bounded synthetic MCMC.
  It is not a saved-results postprocessing command.

Full target-runtime checks are required before a new production run. The Slurm
audit wrapper runs both full checks, using ten one-thread workers, before the
production array. Receipts go to `Results/Engine/` and `Results/Analysis/` here.
The production gate requires matching source hashes, plan hash and runtime;
scheduler dependencies alone do not bypass it.

Historical completed-run source receipts are separate: they are checked against
the original archived bytes in `../Provenance/`, not relabeled as fresh checks
of the reorganized live code. The reference C++ files retain their original
contents. Corresponding R wrappers differ only in locations and explicit
historical/current-receipt handling.
