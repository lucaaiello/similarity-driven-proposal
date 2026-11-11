# Simulation experiments

Each experiment has its own `run.R`, `README.md` and `Results/` folder:

1. `Exact_Gaussian/`
2. `Gaussian_Scaling/`
3. `Correlated_Predictors/`
4. `Exact_Score_Sensitivity/`
5. `Correlated_Score_Sensitivity/`
6. `Dirichlet_Multinomial/Comparison/`
7. `Dirichlet_Multinomial/Fitting_Depth/`

## Shared sampling components

| Location | Components and use |
|---|---|
| Parent folder | `sampler_f_score.cpp`, `sampler_bic_score.cpp`: F/BIC Gaussian samplers shared by the exact-reference and correlated-predictor experiments |
| Parent folder | `gaussian_kernels.cpp`, `gaussian_model.h`: Gaussian comparison/scaling/score-sensitivity kernels and common model calculations |
| `Dirichlet_Multinomial/` | `dm_sampler.cpp`, `dm_numerics.h`, `run_dirichlet_multinomial.R`: shared by the two DM experiments |
| `Validation/` | `gaussian_support.R`, `dm_validation.R`, implementation checks and recorded check results |
| `Figures_and_Tables/` | Shared saved-results summaries, plotting and visual theme |

Commands run from the repository root. `Rscript reproduce/run.R` regenerates
the reported tables/figures without sampling or compilation. Experiment `run.R`
files instead execute MCMC; use `--mode=full` for the manuscript design and
`--mode=quick` for smaller checks. New defaults write to `Results/New_Run/`,
preserving the completed results.

Recorded numerical results are unchanged. Personal paths in saved provenance
metadata have been redacted without changing any nontext values. Source snapshots and embedded
receipts retain their original names/paths as historical provenance. Current
locations are recorded in `Validation/file_locations.csv`. Changed live source
paths require fresh checks before new sampling; old check receipts are not
silently rewritten or presented as current-source validation.
