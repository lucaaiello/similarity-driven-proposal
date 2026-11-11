# Simulation implementation checks

These scripts can compile C++ and execute bounded synthetic chains. They are
not required for saved-results reproduction. Run them explicitly before new
sampling. Recorded original check outputs are under `Results/`; new defaults
use `Results/New_Run/` and do not replace those receipts.

Gaussian order:

```sh
Rscript Simulations/Validation/check_gaussian_kernels.R
Rscript Simulations/Validation/check_gaussian_updates.R
Rscript Simulations/Validation/check_graph_swaps.R
Rscript Simulations/Validation/check_marginal_scores.R
```

Synthetic DM checks:

```sh
Rscript Simulations/Validation/check_dm_numerics.R
Rscript Simulations/Validation/check_dm_adaptation.R
Rscript Simulations/Validation/check_dm_laplace.R
```

Additional coefficient, cache and adaptation checks are provided separately.
Changed source paths require fresh checks before any new MCMC. The original
saved experiment results do not need rerunning to regenerate their tables/figures.
