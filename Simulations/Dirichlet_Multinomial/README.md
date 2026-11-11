# Dirichlet-multinomial simulation experiments

The shared files in this folder support two experiments:

- `Comparison/`: comparison of uniform, fixed/adaptive similarity, graph-guided
  swaps and Gagnon-Laplace. Its completed runs are in `Comparison/Results/`.
- `Fitting_Depth/`: auxiliary fitting budgets one, five and ten, with the same
  coefficient-proposal fitting cap. Full retained fits and parameter diagnostics
  are in `Fitting_Depth/Results/`.

`dm_sampler.cpp` and `dm_numerics.h` implement the common sampler/numerics.
`run_dirichlet_multinomial.R` is the common experiment driver. Each experiment's
`run.R` calls it with the corresponding settings. These are shared DM components;
the Gaussian experiments do not use them. Fit validation is shared separately
under `Simulations/Validation/`.

From the repository root, reproduce the reported outputs without MCMC using
`Rscript reproduce/run.R`. For new sampling, read the experiment README and run
the current-source checks first. Completed results are preserved.
