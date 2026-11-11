FRESH 50K ADAPTIVE TEMPERING: THREE SAMPLERS, FOUR STARTS EACH
===========================================================
This is a new experiment and output tree. Original sources, plans and scientific
results are not overwritten. The verified persistent computational copy is used.
No runtime speedup is assumed from that copy; 10 workers replace four workers.

DESIGN
12 independent ladders, each with 10 replicas and 10 one-thread PSOCK workers.
Methods: online_similarity_flip, online_similarity_local, uniform_flip.
All methods retain the same joint alpha/beta update, warm starts for finite-budget
coefficient optimization, and beta-only mixture. Local swaps within a regression
model are enabled only for online_similarity_local. These predictor moves are
distinct from the temperature exchanges, which are enabled for all three methods.
Coefficient/score fitting caps stay 25/1. Priors and likelihood stay unchanged.

50,000 sweeps PER replica; burn25,000; no thinning; 25,000 cold draws per ladder.
Iterations25,001 through50,000 are retained. Both online adaptations remain active
through iteration50,000. The burn-in cutoff does not freeze the online adaptations.
The local mathematical kernels and finite fitting caps are unchanged.
Ten temperatures are not ten independent posterior replications. Four independent
cold chains per method support the cross-start convergence assessment.

EXCHANGES
Before each set of local sweeps, attempt all nine adjacent edges exactly once,
in a newly sampled random order. Recompute the MH probability after each attempt
using the then-current physical states and the fixed current temperature ladder.
Uniformly mixing every ordering with its reverse preserves the fixed-ladder
joint target. Priors are NOT tempered; only the likelihood is powered.
There are 450,000 exchange attempts per ladder, including50,000 on the cold edge.
Track walker endpoint visits after every exchange attempt.
The MMV log-gap temperature update remains once per iteration, with target0.234.
The existing effective-neighbor/diversity lambda update remains separate/online.
Both adaptations continue through the full run with their existing decreasing gains.

FRESH INITIALIZATION
Null/5%/10%/15% indicators are generated afresh, using exactly round(rate*P*J)
active entries. The same initial indicator configuration is used across methods
for a given start. Every replica in a ladder starts from THAT ladder's start.
All alpha/beta=0, g=1, Sigma_alpha=I; counters/windows, ladder and RNG reset.
No saved posterior endpoint, fitted coefficient or learned covariance is imported.
Global intercept covariance stays I apart from the existing temperature multiplier;
no new intercept-covariance adaptation is introduced. Category joint updates remain.

TASK MAP
1 flip/null;2 local/null;3 uniform/null;
4 flip/5%;5 local/5%;6 uniform/5%;
7 flip/10%;8 local/10%;9 uniform/10%;
10 flip/15%;11 local/15%;12 uniform/15%.
All12 tasks queue together; concurrency cap6, 10 cores/task, up to60 requested
cores concurrently if resources permit. Each task waits for all10 cores on ONE
node (plus40G memory); it does not start partially allocated or use fewer workers.
The first six task IDs include null and5% starts for all3 samplers.
The scheduler chooses dispatch times; this does not enforce simultaneous
starts, one active task per sampler, or equal progress. Uniform may finish earlier.
Given the earlier report of two36-core highp nodes, six10-core ladders can fit
as three per node if cores and memory are available.
Other users and memory/queue constraints may reduce simultaneous starts.
No60-core aggregate reservation is required to begin: eligible tasks start
individually as their10-core allocations become available. h_data4G means40G/ladder.
Each ladder requests120h. Runtime is a request, not a completion prediction.
Extrapolating the previous4-worker5k pilot to10 workers and50k gives roughly83h
for an informed ladder before fresh-state, hardware and scheduling differences.
The120h request is plausible headroom, not a guarantee of completion.

AUDIT AND RECOVERY
The first cluster job runs the existing full computational-equivalence audit
with10 workers, then the new all-edge finite-state/stationarity and bounded
synthetic serial-reference/parallel-candidate/restart checks for ALL3 methods.
Every production task checks the matching successful new audit before sampling,
even if a failed audit job releases its scheduler hold. Summary requires all12
matching final checkpoints, reconstructs retained endpoint targets independently,
and reports fitted disagreements, rank folded split Rhat, PIPs, costs and flow.
No automatic convergence or publishable association claim is made.

Full-state checkpoints every500 sweeps. Keep12.5k/25k/37.5k/50k milestones plus the
latest rolling checkpoint; discard only a superseded rolling checkpoint in THIS
new task output tree after its successor is committed. Old experiments are untouched.
Resume by resubmitting an interrupted task with the same inputs/runtime/plan.
A surviving task_lock requires verifying that its original job is dead first.
The submit script refuses a duplicate submission record. Do not re-run it to monitor.

LOCAL PREPARATION (metadata/packaging only)
powershell.exe -NoProfile -ExecutionPolicy Bypass -File New_simulations\dm_apt50k_fresh_v1\prepare_bundle.ps1

TRANSFER FROM LOCAL POWERSHELL
scp New_simulations\dm_apt50k_fresh_v1_bundle.tgz New_simulations\dm_apt50k_fresh_v1_bundle.tgz.sha256 laiello@dtn.hoffman2.idre.ucla.edu:
ssh laiello@hoffman2.idre.ucla.edu

ON HOFFMAN2: ONE SUBMISSION, AUDITS THEN AUTOMATIC PRODUCTION THEN SUMMARY
(
set -euo pipefail
cd "$HOME"
sha256sum -c dm_apt50k_fresh_v1_bundle.tgz.sha256
work_dir="$HOME/jcgs_dm_apt50k_fresh_v1"
test ! -e "$work_dir" || { echo "STOP: workspace already exists; inspect it."; exit 1; }
mkdir -p "$work_dir"
tar --no-same-permissions -xzf "$HOME/dm_apt50k_fresh_v1_bundle.tgz" -C "$work_dir"
cd "$work_dir"
bash New_simulations/dm_apt50k_fresh_v1/hoffman2/submit.sh
)

MONITOR
cd "$HOME/jcgs_dm_apt50k_fresh_v1"
qstat -u laiello
cat New_simulations/dm_apt50k_fresh_v1/hoffman2/submitted_jobs.txt
Use the recorded ladder job ID in:
tail -n8 New_simulations/dm_apt50k_fresh_v1/hoffman2/logs/ladder.JOB_ID.TASK_ID.log

COMPLETED SUMMARY DOWNLOAD FROM LOCAL POWERSHELL
New-Item -ItemType Directory -Force New_simulations\outputs\dm_apt50k_fresh_v1 | Out-Null
scp -r laiello@dtn.hoffman2.idre.ucla.edu:jcgs_dm_apt50k_fresh_v1/New_simulations/outputs/dm_apt50k_fresh_v1/summary New_simulations\outputs\dm_apt50k_fresh_v1
