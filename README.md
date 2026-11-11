# Similarity-Driven Proposals for MCMC Algorithms on Discrete Spaces

Code, data and saved results accompanying the manuscript.

**You can reproduce the reported tables and figures without rerunning the
samplers.** The saved-results workflow runs locally and needs no compiler or
cluster access.

## Quick start

### 1. Get the files and saved results

If you received `reproduction_package.zip`, extract it. The archive contains
the complete data and saved results.

For a GitHub clone, install [Git LFS](https://git-lfs.com/) before cloning:

```sh
git lfs install
```

Clone using the URL on the repository's GitHub page. Open a terminal in the
downloaded repository folder, then retrieve the large saved result:

```sh
git lfs pull
```

The application file `Real Data Analysis/Results/summary.rds` should be about
232 MiB. A small text file in its place is a Git LFS pointer, not the result.
GitHub's **Download ZIP** may contain pointers; use Git LFS or the complete
reproduction archive in that case.

### 2. Install the R packages

The workflow was tested with **R 4.4.1**. Run this once in R:

```r
install.packages(c("ggplot2", "posterior"), repos = "https://cloud.r-project.org")
```

Allow several GB of free RAM. R must support Cairo PDF output; you can check
this with `capabilities("cairo")` in R.

### 3. Reproduce the tables and figures

From the **repository root**—the folder containing this README and the four
folders listed below—run:

```sh
Rscript reproduce/run.R
```

On Windows, if `Rscript` is not recognized, use your R installation's full
path. For example, in PowerShell:

```powershell
& 'C:\Program Files\R\R-4.4.1\bin\Rscript.exe' reproduce/run.R
```

Adjust the installation path if your R version or location differs. The
repository itself can be stored anywhere; R scripts use paths relative to
the repository root.

The default run usually takes one to a few minutes on a laptop. It creates
`reproduced/` and leaves the supplied results untouched.

## Find the outputs

| Location inside `reproduced/` | Contents |
|---|---|
| `figures/` | Main and supplementary figures in PDF and PNG formats; the model diagram is PDF only |
| `tables/` | Numerical table summaries as CSVs, plus formatted table exports |
| `table_verification.csv` | Agreement checks for all 16 numerical tables |
| `status.txt` | `COMPLETE` when every step has passed |
| `logs/`, `timings.csv` | Step logs and runtimes |
| `sessionInfo.txt`, `package_versions.csv` | Software environment used |

The exact figure and table mapping is in
[docs/REPRODUCTION.md](docs/REPRODUCTION.md).

To repeat the workflow, choose a new output folder:

```sh
Rscript reproduce/run.R --output=reproduced_second_check
```

Existing output folders are preserved. If a step fails, inspect its log in
`logs/`; `status.txt` remains `INCOMPLETE`.

## Repository guide

| Folder | What you will find |
|---|---|
| [Simulations/](Simulations/README.md) | Simulation experiments, their samplers and saved results. Each experiment has its own instructions and `Results/` folder. |
| [Real Data Analysis/](Real%20Data%20Analysis/README.md) | Diet–microbiome data, the application sampler, saved results, diagnostics and optional cluster scripts. |
| [reproduce/](reproduce/) | The shared workflow for rebuilding and checking all reported tables and figures. |
| [docs/](docs/) | Data descriptions, detailed reproduction instructions and instructions for new sampling runs. |

## More options

The default workflow uses saved parameter diagnostics. To recompute them from
the saved fits too, use:

```sh
Rscript reproduce/run.R --output=reproduced_full_diagnostics --recompute-parameter-diagnostics=true
```

This takes additional time but still runs no sampler or optimizer.

For new simulation or application runs, follow
[docs/RERUN_EXPERIMENTS.md](docs/RERUN_EXPERIMENTS.md). Those commands execute
MCMC and may require compilation; they are optional for reproducing the paper.

To build a standalone ZIP containing the code, data and saved results, use
Python 3.9 or later:

```sh
python reproduce/prepare_archive.py
```

The verified archive and its checksum are written to `dist/`. Existing archives
are preserved; the script does not upload anything.

## Further information

- [Data and saved-result descriptions](docs/DATA.md)
- [Figure/table mapping and calculation details](docs/REPRODUCTION.md)
- [Instructions for running new experiments](docs/RERUN_EXPERIMENTS.md)
- [Tested software environment](docs/validated_environment.txt)
- [Licensing status](LICENSING.md) — a software license has not yet been selected;
  third-party data retain their own terms.
