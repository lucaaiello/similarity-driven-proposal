# Figure and table map

Run `Rscript reproduce/run.R` from the repository root. No sampling, optimization,
compilation or remote computation is involved.

## Figures

All filenames below are under `reproduced/figures/`, with PDF and PNG versions
unless stated otherwise.

| Figure | Output basenames |
|---|---|
| Main 1 | `simulation_exact_reference_accuracy_ggplot`, `simulation_exact_reference_efficiency_ggplot` |
| Main 2 | `simulation_correlated_local_moves_movement_ggplot`, `simulation_correlated_local_moves_efficiency_ggplot` |
| Main 3 | `simulation_auxiliary_score_alignment_ggplot`, `simulation_auxiliary_score_tuning_ggplot`, `simulation_auxiliary_score_pips_ggplot` |
| Main 4 | `simulation_dm_comparison_recovery_ggplot`, `simulation_dm_comparison_efficiency_ggplot` |
| Main 5 | `real_data_spreads` |
| S1 | `simulation_highdim_scaling_throughput_ggplot`, `simulation_highdim_scaling_efficiency_ggplot` |
| S2 | `DM_diagram_2.pdf` (static diagram copied from the source assets) |
| S3 | `real_data_pip_stability` |
| S4 | `real_data_associations` |
| S5 | `real_data_retained_parameters` |

Synthetic figures use `make_manuscript_simulation_figures_ggplot.R`; application
figures use `make_manuscript_real_data_figures.R` and its ggplot helper. All share
`manuscript_plot_theme.R`. The additional combined plots in `supporting_figures/`
are base-graphics exports supporting the numerical-summary workflow, not the
manuscript's current panel layout.

## Tables

Files are under `reproduced/tables/`. `reproduce/tables.R` calculates every
computational summary; `check_tables.R` compares its entries with the current
explicit TeX tables at their reported precision.

| Table | Output |
|---|---|
| Main 1, DM transition definitions | `tab_DM_methods.tex` (static) |
| Main 2, synthetic DM diagnostics | `main2_dm_computation.csv` |
| Main 3, application comparison | `main3_application_comparison.csv` |
| Main 4, supported associations | `main4_application_associations.csv` |
| S1, designs | `tab_simulation_designs.tex` (static) |
| S2, proposal definitions | `tab_simulation_methods.tex` (static) |
| S3, exact Gaussian results | `S3_exact_gaussian.csv` |
| S4, scaling | `S4_scaling.csv` |
| S5, scaling diagnostics | `S5_scaling_diagnostics.csv` |
| S6, correlated-model movements | `S6_correlated_moves.csv` |
| S7, exact score sensitivity | `S7_exact_score.csv` |
| S8, correlated score sensitivity | `S8_correlated_score.csv` |
| S9, DM dataset diagnostics | `S9_dm_dataset_diagnostics.csv` |
| S10, DM parameter diagnostics | `S10_parameter_diagnostics.csv` |
| S11, fitting-depth sensitivity | `S11_fitting_depth.csv` |
| S12, DM recovery | `S12_dm_recovery.csv` |
| S13, DM computation | `S13_dm_computation.csv` |
| S14, application chains | `S14_application_chains.csv` |
| S15, application scalar diagnostics | `S15_application_scalars.csv` |

Static tables describe designs/methods and are retained as table-only reference
snippets in `reproduce/reference/`; they are not estimates. That directory also
contains the model diagram. Full manuscript documents are not distributed here.

## Aggregation and presentation

- Average chains within datasets/scenarios before means and Monte Carlo standard
  errors across datasets. Chains are not treated as independent data replicates.
  Exact Gaussian summaries cover 180 scenario/dataset combinations.
- Robust ESS/s is the smaller of the two monitored scalar ESS/s values, calculated
  before averaging. It is not the minimum of two final aggregate means.
- CSVs retain unrounded proportions. S3/S7 posterior errors use percentage points
  in print; S12 false-positive rates use percentages.
- PIP ranges first compare chains for each predictor. Application Main 3 takes
  each genus's maximum spread and then its median across genera.
- Default S10 reproduction reaggregates saved parameter-level diagnostics.
  `--recompute-parameter-diagnostics=true` recalculates them from full fits.
  Excluded coefficients remain zero; constant-chain diagnostics are flagged, not
  treated as converged. Threshold counts concern available finite diagnostics.
- S5 guided-method acceptance means/MCSEs were rounded to three decimals before
  printing a trailing zero. S7 final-lambda summaries were rounded to three
  decimals before their two-decimal display. The verifier explicitly preserves
  these original intermediate-rounding conventions; CSVs remain full precision.
- Main 4 uses empirical PIP >= 0.5 in every one of the 12 cold chains. This shared
  support criterion is not a Bayesian FDR rule.
- Numerical values are calculated from saved results, not copied from TeX. The
  verifier reads printed numbers only for the final comparison. Tables remain
  explicitly written in the manuscript sources, without generated inputs.

## Runtime, commands and verification

The default typically takes one to a few minutes on a recent laptop. Measured
validation timings are in `validation_timings.csv`. Allow several GB of free RAM
for the 232 MiB application summary and plotting/diagnostic arrays. Full
parameter-diagnostic recomputation takes additional time; no sampler is run.

Each run writes logs, versions, timings, `table_verification.csv` and a completion
status. Failures leave `INCOMPLETE`; inspect the relevant log and use a new
destination after fixing the cause. The master command refuses overwriting.

Individual figure-only commands with separate destinations:

```sh
Rscript New_simulations/make_manuscript_simulation_figures_ggplot.R --image-dir=reproduced_figures
Rscript New_simulations/make_manuscript_real_data_figures.R --output-root=reproduced_application
```

Their legacy defaults write to `JCGS_submission/Images`, so specify destinations
when preserving those assets. Additional assessment code is retained in
`dm_apt50k_postprocessing/analyze_summary.R`; it is not needed for the displays.

Input MD5s are in `reproduce/inputs_md5.csv`; original experiment receipts are
also retained. `.gitattributes` prevents checkout line-ending conversion from
invalidating source hashes. The journal archive adds per-member SHA-256 hashes.
Numerical reproduction does not imply seeded trajectory identity across compilers
or new evidence of convergence. Historical sampling times are observations,
not runtimes of the local postprocessing command.
