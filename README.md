# Similarity-Driven Proposals for MCMC Algorithms on Discrete Spaces

Reproducibility materials for the manuscript by Luca Aiello, Raffaele Argiento,
Alexandros Beskos and Maria De Iorio.

**The saved results are included.** Rebuilding the reported tables and figures
requires no MCMC, C++ compilation, Hoffman2 access or cluster submission.

## Get the complete materials

For a Git clone, install [Git LFS](https://git-lfs.com/) first:

```sh
git lfs install
git clone https://github.com/lucaaiello/similarity-driven-proposal.git
cd similarity-driven-proposal
git lfs pull
```

The 232 MiB application `summary.rds` is stored through Git LFS.
A small text pointer is not the saved result. Alternatively, extract the complete
`reproduction_package.zip` supplied with the article: it contains the actual
result objects and needs no Git LFS. A GitHub source ZIP may contain LFS pointers.

## Software

Saved-results reproduction was tested with **R 4.4.1**, `ggplot2` and `posterior`.
No compiler is required. Install the packages in R:

```r
install.packages(c("ggplot2", "posterior"), repos = "https://cloud.r-project.org")
```

The tested versions/environment are in [docs/validated_environment.txt](docs/validated_environment.txt).
Each reproduction run also records its own environment. Cairo PDF support is
required (`capabilities("cairo")`); the tested Windows R installation provides it.
PDF metadata/font rendering may differ across systems even when numbers agree.

## Reproduce every table and figure

Open a terminal **in the repository root**, then run:

```sh
Rscript reproduce/run.R
```

On Windows, if `Rscript` is not on your PATH:

```powershell
& 'C:\Program Files\R\R-4.4.1\bin\Rscript.exe' reproduce/run.R
```

This creates a new `reproduced/` directory and refuses to overwrite one. To repeat
the check, use `Rscript reproduce/run.R --output=reproduced_second_check`.

| Output | Contents |
|---|---|
| `status.txt` | `COMPLETE` only after every step succeeds |
| `figures/` | Main/supplementary PDF and PNG panels; model diagram PDF |
| `tables/` | Full-precision CSVs for all computational tables; static definition tables as TeX |
| `tables/application_latex/` | Formatted application table exports |
| `table_verification.csv` | Checks of all 16 computational tables against the manuscript |
| `timings.csv`, `logs/` | Step runtimes, exit status and logs |
| `sessionInfo.txt`, `package_versions.csv` | Environment used for this run |

The workflow verifies saved-input checksums, recalculates the numerical summaries
and compares every computational table with its explicitly typeset counterpart.
See [docs/REPRODUCTION.md](docs/REPRODUCTION.md) for the exact figure/table map,
aggregation and rounding conventions. Original results and manuscript assets are
not overwritten. The manuscript tables remain explicitly written in TeX; no
generated table is automatically inserted into them.

The default uses saved parameter-level diagnostics for Table S10. To recompute
those diagnostics from saved full fits as well:

```sh
Rscript reproduce/run.R --output=reproduced_full_diagnostics --recompute-parameter-diagnostics=true
```

This also runs no sampler or optimizer, but takes more time. The default workflow
takes roughly one to a few minutes on a laptop and requires several GB of free
RAM; measured timings are in [docs/validation_timings.csv](docs/validation_timings.csv).

## Organization and saved results

| Location | Purpose |
|---|---|
| `reproduce/` | Master command, table calculations/checks, input manifest, archive builder |
| `docs/` | Result map, data provenance, environment and optional rerun instructions |
| `JCGS_submission/` | Current main/supplement sources, PDFs, bibliography and figure assets |
| `New_simulations/` | Experiment drivers, helpers, audits and saved outputs |
| `Simulations/` | Gaussian C++ samplers used by the final experiments |
| `Real Data Analysis/Data/` | Supplied counts, scaled dietary covariates and labels |

Saved-result directories below are under `New_simulations/outputs/`:

| Analysis | Directory |
|---|---|
| All five matched Gaussian experiments | `matched_gaussian_20261007_163020/` |
| Synthetic DM uniform/similarity methods | `dm_online_adaptive_primary_20datasets/` |
| Synthetic DM Gagnon--Laplace comparator | `dm_gagnon_laplace_20datasets/` |
| Combined synthetic DM comparison | `dm_gagnon_comparison/` |
| DM fitting-depth study, including full fits | `dm_fitting_depth_20261008_104824/` |
| Final application: 12 completed ladders | `dm_apt50k_fresh_slurm_v1/summary/` |

The application summary contains 25,000 retained cold draws per ladder after
25,000 burn-in sweeps, with no thinning. There are four independent cold chains
per method. The ten temperature replicas are not ten independent posterior
chains. All 30 genera enter the model; summaries without Prevotella use those
same full-model draws and do not involve a refit. Remaining diagnostic
qualifications are preserved. Reproduction checks calculations, not new
convergence. Full remote continuation checkpoints are not all included locally.

Paths are preserved because original provenance checks depend on them.
`dm_apt_runtime_optimization` and retained `dm_adaptive_tempering` modules are
dependencies of the final `dm_apt50k_fresh_v1` workflow, not alternative pilot
instructions. See [docs/DATA.md](docs/DATA.md) and
[docs/RERUN_EXPERIMENTS.md](docs/RERUN_EXPERIMENTS.md). Sampling is optional and
much more expensive than the saved-results workflow.

## Build the documents (optional)

With a TeX installation containing `latexmk`:

```sh
cd JCGS_submission
latexmk -pdf -interaction=nonstopmode -halt-on-error main_r1.tex supplements_r1.tex
```

The bibliography style and cross-reference configuration are included. The
validated main/supplement PDFs have 35/44 pages. Building them runs no analyses.

## Journal archive and release

Using Python 3.9 or later:

```sh
python reproduce/prepare_archive.py
```

This builds `dist/reproduction_package.zip` and its SHA-256 checksum, verifies
every member, and excludes Git history, local backups, regenerated results,
response/cover letters and internal housekeeping. It rejects LFS pointers.
It performs no upload, commit or publication.

Licensing is intentionally undecided: [LICENSING.md](LICENSING.md). Third-party
data retain their own terms. Before release, select the code license with the
coauthors, confirm data redistribution terms, publish the prepared revision,
and archive a fixed release with a persistent identifier. Supply the complete
code/data archive to the journal as well as the repository link.

See [docs/JOURNAL_CHECKLIST.md](docs/JOURNAL_CHECKLIST.md), the
[JCGS instructions](https://www.tandfonline.com/action/authorSubmission?show=instructions&journalCode=ucgs20)
and the [ASA discussion of JCGS reproducibility](https://magazine.amstat.org/blog/2017/10/01/dianne-cook/).
The current submission instructions should be checked at upload; this package
does not claim journal certification.
