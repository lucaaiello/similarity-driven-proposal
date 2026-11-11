# Application figures and tables

`create_figures.R` reads the completed `Results/summary.rds` and prepared inputs,
verifies the historical source receipts against the original source archive,
and regenerates numerical summaries, formatted table exports and figures.
`plot_figures.R` draws the panels using the same theme as the simulations.

From the repository root:

```sh
Rscript "Real Data Analysis/Figures_and_Tables/create_figures.R" --output-root=reproduced_application
```

The full checked workflow is `Rscript reproduce/run.R`. Neither command advances
a chain, performs fitting or compiles C++. Completed saved results remain intact.
The optional `--output-root` controls all generated destinations; its default
is `reproduced_application`, not the private submission folder.
