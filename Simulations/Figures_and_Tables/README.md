# Shared simulation figures and tables

- `plot_results.R`: manuscript-quality simulation figure panels, from saved results.
- `summarize_results.R`: chain/dataset aggregations and supporting figure exports.
- `plot_theme.R`: shared visual theme, also used for the real-data figures.

These scripts perform postprocessing, not MCMC. Run from the repository root.
The full checked workflow is `Rscript reproduce/run.R`. To generate only the
simulation figures in a separate directory:

```sh
Rscript Simulations/Figures_and_Tables/plot_results.R --image-dir=reproduced_simulation_figures
```
