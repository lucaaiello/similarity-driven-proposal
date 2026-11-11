# Optional sampling runs

These are **not needed** for `Rscript reproduce/run.R`. They compile C++, execute
MCMC or submit jobs only when explicitly invoked. Always use fresh destinations.

## Gaussian experiments

Install Rtools matching R 4.4 on Windows and the packages checked by the drivers,
including Rcpp, RcppArmadillo and RcppEigen. From the repository root:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File New_simulations/run_matched_gaussian_revision.ps1
```

This runs validation gates and the five matched Gaussian experiments. See
`New_simulations/MATCHED_GAUSSIAN_REVISION.md` for settings/options. Its historical
runtime estimates are approximate. Seeds are set in individual drivers and
recorded in provenance. BLAS/OpenMP threads are limited to one for comparisons.

## Synthetic DM

Compiler/package requirements are in the audits and
`New_simulations/hoffman2/install_pinned_packages.R`. The principal design is
20 datasets, three chains, 5,000 sweeps, burn 2,500, thinning five. Examples:

```sh
Rscript New_simulations/run_dm_numerical_stability_audit.R
Rscript New_simulations/run_dm_online_adaptation_audit.R
Rscript New_simulations/run_dm_gagnon_laplace_audit.R
Rscript New_simulations/run_dm_online_adaptive_ablation.R --mode=full --stage=core --datasets=20 --chains=3 --runs=5000 --burn=2500 --thin=5 --methods=uniform_flip,optimized_similarity_fixed,online_similarity_flip,online_similarity_local --output=New_simulations/outputs/new_dm_primary
Rscript New_simulations/run_dm_online_adaptive_ablation.R --mode=full --stage=core --datasets=20 --chains=3 --runs=5000 --burn=2500 --thin=5 --methods=gagnon_laplace_rj --output=New_simulations/outputs/new_dm_gagnon
```

The fitting-depth study uses datasets 1–5, fixed similarity, score fitting
budgets 1/5/10, `--save-full-fits=true`, `--cache-score-fits=true`, three chains,
5,000 sweeps, burn 2,500 and thinning one. Exact settings are retained in each
depth's manifest. Run the required audits first. Sampling is much more expensive
than postprocessing and is not automatically invoked by the reproduction script.

## Application production

The final entry point is `New_simulations/dm_apt50k_fresh_v1/slurm/submit.sh` in a
new Hoffman2 workspace, after installing the pinned environment and passing the
matching audits. It requests 12 ladders, ten one-thread workers each, at most six
concurrent ladders, and 120 hours per task. Read the Slurm README/environment
scripts first. Original account/library settings are retained; adapt scheduler
wrappers for another environment, update their manifests and rerun the gates.
Do not bypass scientific integrity checks. Legacy Grid Engine wrappers remain
only because the frozen bundle references them.

Design: 50,000 sweeps per replica, 25,000 burn-in, no thinning; all adjacent
exchange attempts per sweep; fresh null/5%/10%/15% starts with continuous states
and tuning reset. Likelihoods are tempered; priors are not. Earlier informed
ladders required about two days each, with queue time additional. Hardware and
allocation strongly affect runtime.

Local experiments used a Dell XPS 17 9700, Intel Core i7-10750H (six physical
cores, 12 logical processors), approximately 16 GB RAM. The final application
ran on UCLA Hoffman2 with R 4.4.1, Rcpp 1.0.12, RcppArmadillo 15.2.3.1 and
RcppEnsmallen 0.2.22.1.1. Worker BLAS/OpenMP threads were limited to one. Exact
trajectories/timings can depend on compiler, numerical libraries and hardware.
