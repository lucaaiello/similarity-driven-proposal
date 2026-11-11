# Synthetic Dirichlet–multinomial comparison

Uniform, fixed/adaptive similarity, graph-guided swaps and Gagnon–Laplace. Main Figure 4, Table 2; Supplementary Tables S9, S12–S13.

- `run.R`: experiment entry point. This executes MCMC; it is not needed to reproduce saved results.
- `Results/`: completed experiment results. Files are preserved from the original runs.
- Shared DM samplers, numerics and the experiment driver are in the parent `Dirichlet_Multinomial/` folder.

From the repository root, regenerate all reported tables/figures without MCMC:

```sh
Rscript reproduce/run.R
```

To rerun this experiment, first run the required checks in `Simulations/Validation`, then invoke:

```sh
Rscript Simulations/Dirichlet_Multinomial/Comparison/run.R --mode=full
```

New sampling outputs go to `Results/New_Run`; choose another `--output` for another run.
Completed results are never overwritten by these defaults. Numerical receipts and
source snapshots within `Results` retain their original recorded paths/names;
these describe the completed run, rather than the current directory structure.
