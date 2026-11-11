# Auxiliary-score sensitivity with correlated predictors

Comparison of conditional and marginal guidance in the correlated design. Main Figure 3 (last panel); Supplementary Table S8.

- `run.R`: experiment entry point. This executes MCMC; it is not needed to reproduce saved results.
- `Results/`: completed experiment results. Files are preserved from the original runs.
- Shared samplers and numerical helpers are in the parent `Simulations/` folder.

From the repository root, regenerate all reported tables/figures without MCMC:

```sh
Rscript reproduce/run.R
```

To rerun this experiment, first run the required checks in `Simulations/Validation`, then invoke:

```sh
Rscript Simulations/Correlated_Score_Sensitivity/run.R --mode=full
```

New sampling outputs go to `Results/New_Run`; choose another `--output` for another run.
Completed results are never overwritten by these defaults. Numerical receipts and
source snapshots within `Results` retain their original recorded paths/names;
these describe the completed run, rather than the current directory structure.
