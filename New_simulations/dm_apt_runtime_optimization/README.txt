PERSISTENT TEMPERING COMPUTATIONAL COPY AND REGRESSION
=====================================================
This phase changes computation only. The existing ONE uniformly random adjacent
exchange per iteration,10-replica design,0.234 ladder target,diversity adaptation,
proposal settings and posterior target are unchanged. More-frequent exchanges
and the fresh-start long production experiment are NOT implemented or launched
by this verification workflow. Existing sources/results/checkpoints are immutable.

NATIVE OPTIMIZATION
DM_persistent.cpp retains immutable X/Y,validated neighborhood,row totals and
uniform probabilities in a context. It reuses physical-state views and replaces
the full production function output with one compact step result. It removes
the redundant zero-state initialization,second graph/data setup,production
retained buffers,large config/continuation construction and banner/progress work.
The initial physical rebase and final full likelihood refresh are preserved.
Optimizer caps,finite-iteration starts,exact cache keying and random-draw order
are unchanged. Prefix numerical helpers and live sweep are derived from a
byte-identical reference snapshot;native_generation.rds records the provenance.

Native contexts are reconstructible workspaces,not serialized Markov state.
Physical states,tuning,temperature/diversity windows,and RNG stay in the same
R ladder state. Invalid/deserialized native pointers are rejected. Every step
accepts the entire physical state,so after a swap no stale replica cache is used.
The first parity implementation deliberately retains apt_score unchanged in R.

SUPPORTED SCOPE
The3 existing similarity/uniform methods under current production priors10/10,
a1/b9 and the fixed within-beta mixture plus joint scale2 are supported.
This compact interface is not a replacement for arbitrary sampler configurations.
Fresh null/5/10/15% synthetic fixtures reset alpha/beta=0,g=1,Sigma_alpha=I.
Saved real-data states are bounded TEST FIXTURES only,not scientific long-run
initialization. No convergence,posterior-accuracy or new start-policy claim.
Future fresh production needs its explicit warmup policy;we do not silently
add covariance adaptation during this computational-parity phase.

VERIFICATION
parse:parse sources only.
compile:compile reference/new copies;no optimizer/MCMC.
fixed:construct contexts and compare fixed-state helpers;no optimizer/MCMC.
quick:bounded synthetic worker/whole-ladder MCMC,serialized restarts and worker
  scheduling checks.24sweeps,3replicas for the ladder case.
full:80sweeps,10replicas for synthetic whole-ladder comparisons,and6 paired real
  single-step fixtures. Not a real-data chain or production experiment.
Every boundary compares adaptation,physical/coefficient states,temperatures,
swap decisions,labels,diagnostics and RNG. Indicators/RNG/counters must match
exactly;numerical fields use tolerance1e-10. Per-step diagnostic projection
counters are tested,not hidden. Restarts rebuild native pointers while preserving
serialized R state. Parallel scheduling is checked for selected worker count.
Protected source/checkpoint/summary hashes are checked before and after.

LOCAL COMMAND
powershell.exe -NoProfile -ExecutionPolicy Bypass -File New_simulations\dm_apt_runtime_optimization\run_verification.ps1 -Mode full -Workers 2
Run only when bounded verification is intended. It never starts production or
submits jobs. Results go to a new dm_apt_optimization_* output directory.

OPTIONAL HOFFMAN2 CHECK
prepare_bundle.ps1 creates a verification-only archive with the immutable reference
source graph and small diagnostic inputs. It invokes no computation or submission.
Extract it in a NEW jcgs_apt_runtime_optimization_v1 workspace, create
New_simulations/dm_apt_runtime_optimization/hoffman2/logs, then submit verify.qsub
from that workspace. That single optional job tests 10 one-thread workers in the
cluster environment. No scientific sampling array is included. No cluster job
has been submitted by this implementation. Do not overwrite old workspaces.

PERFORMANCE INTERPRETATION
Compilation is excluded from timing probes;engine order rotates. Bounded real
single-step timings are descriptive,not a whole-chain cluster speedup guarantee.
The expensive neighborhood/finite-budget coefficient fits remain. Any claimed
whole-run speedup needs a matched timing experiment after parity checks pass.
Changing4workers to5or10 can reduce wall time without reducing total CPU work;
its seeded scheduling parity is covered by the selected worker-count audit.

LOCAL RESULT, 2026-10-29
Full verification passed: 148 checks, 10-replica bounded synthetic ladders,
2 PSOCK workers, exact indicators/RNG/swap decisions, numerical tolerance 1e-10,
and six paired full-size real-data single-step fixtures. Serialized restarts
reconstructed the native contexts. Original sources and saved results retained
their hashes. See outputs/dm_apt_optimization_full_20261029_121219_26796/report.txt.

These fixtures did NOT demonstrate a useful speedup: similarity local steps
were about 4 seconds with either engine, and similarity flip steps about
3.4-3.7 seconds. This copy has passed equivalence testing; it is not evidence
of faster production. Avoid extrapolating a shorter long-run runtime from it.
