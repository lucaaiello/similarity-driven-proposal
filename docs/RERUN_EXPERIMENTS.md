# Optional sampling runs

These commands execute MCMC and can compile C++. They are not needed for
`Rscript reproduce/run.R`, which uses the saved results. Run from the repository
root with Rcpp, RcppArmadillo, RcppEigen and the packages checked by each driver;
on Windows use matching Rtools. New results go to each experiment's
`Results/New_Run/` and do not replace the completed results.

## Gaussian experiments

Run the implementation checks first:

```sh
Rscript Simulations/Validation/check_gaussian_kernels.R
Rscript Simulations/Validation/check_gaussian_updates.R
Rscript Simulations/Validation/check_graph_swaps.R
Rscript Simulations/Validation/check_marginal_scores.R
```

Then run the required experiment(s), sequentially for timing comparisons:

```sh
Rscript Simulations/Exact_Gaussian/run.R --mode=full
Rscript Simulations/Gaussian_Scaling/run.R --mode=full
Rscript Simulations/Correlated_Predictors/run.R --mode=full
Rscript Simulations/Exact_Score_Sensitivity/run.R --mode=full
Rscript Simulations/Correlated_Score_Sensitivity/run.R --mode=full
```

Limit OMP/OpenBLAS/MKL threads to one when replicating the original computational
comparisons. Seeds and settings remain in the experiment drivers and retained
manifests. `--mode=quick` selects smaller runs; they are not manuscript results.
Compilation, enumeration and sampling can take hours. Do not rerun MCMC merely
to regenerate the reported figures/tables.

## Synthetic Dirichlet-multinomial experiments

```sh
Rscript Simulations/Validation/check_dm_numerics.R
Rscript Simulations/Validation/check_dm_adaptation.R
Rscript Simulations/Validation/check_dm_laplace.R
Rscript Simulations/Dirichlet_Multinomial/Comparison/run.R --mode=full
Rscript Simulations/Dirichlet_Multinomial/Fitting_Depth/run.R --mode=full
```

The main comparison uses 20 datasets, three chains, 5,000 sweeps, burn 2,500 and
thinning five. Its entry point runs similarity methods and Gagnon-Laplace into
separate result subfolders. The fitting-depth experiment uses datasets 1-5,
three chains, 5,000 sweeps, burn 2,500, no thinning, and auxiliary fitting caps
one/five/ten; the coefficient fitting cap stays three. Full fits are retained.
Both call the shared `Simulations/Dirichlet_Multinomial/run_dirichlet_multinomial.R` driver.

For a newly sampled comparison, call `summarize_comparison.R` with explicit
`--baseline`, `--gagnon` and a fresh `--output` directory. Default postprocessing
reads the original completed experiment. Use new destinations for reruns.

## Application production

The optional entry point is `Real Data Analysis/Cluster/submit.sh` in a new
Hoffman2 workspace. Read the application and Cluster READMEs first. Prepare a
source/input bundle locally with:

```sh
python "Real Data Analysis/Cluster/prepare_bundle.py"
```

This prepares files only. On the cluster, after extracting into a new workspace
and installing the pinned environment, explicitly invoke:

```sh
bash "Real Data Analysis/Cluster/submit.sh"
```

The Slurm chain requires full engine and analysis audits before 12 ladders
(three methods x four starts), ten one-thread workers each, with at most six
concurrent ladders. The retained requests are 40 GB and 120 hours per ladder.
They are recorded settings; adapt the account/module/library paths to the
target site. New states and summary go to `Results/New_Run`, preserving the
completed `Results/summary.rds`. No cluster jobs are invoked by reproduction.

The scientific settings remain 50,000 sweeps per replica, burn-in 25,000, no
thinning, and all adjacent exchange attempts per sweep. Fresh null/5%/10%/15%
starts reset continuous states and tuning. Only the likelihood is tempered.
Completed similarity-based ladders required about two days each; queue time is additional.
The frozen plan and prepared inputs are unchanged. Original source receipts
are verified against the archived original source files; new sampling requires
fresh checks of the current executable paths and runtime. Full old remote
checkpoints are not all included, and are not relabeled as compatible resumes.

Local experiments used a Dell XPS 17 9700, Intel Core i7-10750H (six physical
cores, 12 logical processors), approximately 16 GB RAM. The completed application
used UCLA Hoffman2, R 4.4.1, Rcpp 1.0.12, RcppArmadillo 15.2.3.1 and RcppEnsmallen
0.2.22.1.1, with worker BLAS/OpenMP threads limited to one. Exact trajectories
and timings can depend on compiler, numerical libraries and hardware.
