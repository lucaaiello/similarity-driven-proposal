# Corrected Gaussian reruns

From the repository root, run the following once in PowerShell:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\New_simulations\run_matched_gaussian_revision.ps1
```

The launcher runs four audits, five small pilot experiments, then the five full
Gaussian experiments. It stops on any nonzero R exit code. To run only the audits
and pilots, append `-PilotOnly`. To choose the new results directory, append
`-OutputRoot 'C:\path\to\a\NEW\folder'`. Existing directories are refused.

The launcher fixes OMP/OpenBLAS/MKL thread settings to one for the run and restores
the previous settings afterwards. Keep this machine otherwise idle for the
wall-clock comparisons.

## Changes

- Production Gaussian model acceptance uses the collapsed Beta--binomial
  configuration prior. Conditional phi draws are retained as auxiliary output
  for compatibility, but do not enter the acceptance ratio or log-target trace.
- Fixed F/BIC proposals and their adaptive counterparts both use balanced
  birth/death directions, a direction-preserving 2% uniform mixture, the same
  rank-one engine and proposal caching. Their initial lambda is 0.7; only the
  adaptation flag differs within each comparison.
- The auxiliary scores no longer impose an n-2 cutoff on posterior support.
  Singular/dense OLS fits use a pseudoinverse fallback. F scores are neutral
  when residual degrees of freedom are unavailable; all model configurations
  retain a proper conjugate posterior target.
- Timed comparator runs now include Gram-matrix setup, matching the optimized
  similarity timer. Both timed engines require at least 50 post-burn samples
  before a time-budget stop. Actual elapsed times, not nominal budgets, are used.
- The local-move benchmark retains its original block-change diagnostic and
  additionally records direct substitutions between distinct single-predictor
  block representatives, excluding transitions to/from the zero code.
- Every corrected run records source hashes, the invocation (including settings),
  session information and thread settings. Existing experiment-specific
  manifests/data and chain seeds are preserved.

## Audit scope

The new audit enumerates all 32 states of regular, dense and duplicate-column
five-predictor designs. It reconstructs frozen one-flip MH transition matrices
from C++ proposals/targets and checks them against an independently coded R
posterior. It checks normalization, positive neighbor support, balanced
directions, detailed balance and stationarity. Checks of actual production runs
cover fixed/frozen-adaptive equality, active adaptation, recorded log targets,
proposal-cache refreshes and independence from initial phi. These frozen-kernel
checks do not establish a general online-adaptation ergodicity theorem.

The existing rank-one, swap and marginal-score audits also rerun, with current
source-hash checks before the associated experiments.

## Full design and timing

| Stage | Design | Approximate duration |
|---|---|---|
| Audits + pilots | Four audits, five short runs | 10--20 min |
| Exact Gaussian | Nine scenarios, 20 datasets each, 3 chains, 7 methods; 20,000 iterations / 5,000 burn | 20--35 min |
| Scaling | Five dimensions, 20 datasets, 3 chains, 7 methods; 5 seconds per chain | 3--3.3 h |
| Correlated local moves | 30 datasets, 4 chains, 6 methods; 50,000 iterations / 10,000 burn | 2--2.7 h |
| Exact score sensitivity | 20 datasets, 4 chains; 3 adaptive scores + 15 fixed settings; 1 second per chain | 25--35 min |
| Correlated score sensitivity | 10 datasets, 4 chains, 3 scores; 5 seconds per chain | 10--15 min |
| Figures and tables | Corrected Gaussian data only | Under 2 min |

Allow roughly 6--7 hours sequentially; timings depend on machine load. The old
saved chain-time totals were 15.9 minutes for exact sampling, 2.92 hours for
scaling, and 2.21 hours for local moves. Those totals exclude compilation,
enumeration, data generation and output serialization. The launcher records
actual stage durations in `pipeline_timings.csv`.

## Outputs and manuscript handoff

Each launch creates `New_simulations/outputs/matched_gaussian_YYYYMMDD_HHMMSS`.
The full experiment folders are `exact`, `scaling`, `swap`, `exact_score` and
`correlated_score`. Stage logs and audit/pilot folders sit beside them.

After successful full runs, corrected Gaussian figures are written to
`<run root>/Images` and dataset-level table summaries to
`<run root>/manuscript_simulation_summaries`. The launcher deliberately does not
overwrite `JCGS_submission/Images` or any LaTeX file. Inspect the corrected
results before replacing manuscript assets and reported numbers.

Dirichlet--Multinomial experiments are not rerun by this launcher. Historical
Gaussian experimental variants are not manuscript-ready replacements for these
matched comparisons. The separate lambda-grid scripts now default to balanced
directions and a 2% mixture, but their pilot-selection experiments are not needed
for the five current manuscript reruns.
