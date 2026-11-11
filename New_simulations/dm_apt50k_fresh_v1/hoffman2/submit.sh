#!/bin/bash
set -euo pipefail
root=New_simulations/dm_apt50k_fresh_v1
test -f "$root/inputs/plan.rds"
sha256sum -c "$root/inputs/bundle_files.sha256"
mkdir -p "$root/hoffman2/logs"
record="$root/hoffman2/submitted_jobs.txt"
( set -o noclobber; printf 'FRESH50K: 3 SAMPLERS x4 STARTS x10 REPLICAS\n' > "$record" ) || {
    echo "Submission record exists. Do not submit duplicates."; exit 1;
}
audit=$(qsub -terse "$root/hoffman2/audit.qsub"); aid=${audit%%.*}
printf 'audit %s\n' "$aid" >> "$record"
ladders=$(qsub -terse -hold_jid "$aid" "$root/hoffman2/ladders.qsub"); lid=${ladders%%.*}
printf 'ladders %s\n' "$lid" >> "$record"
summary=$(qsub -terse -hold_jid "$lid" "$root/hoffman2/summary.qsub"); sid=${summary%%.*}
printf 'summary %s\n' "$sid" >> "$record"
cat "$record"
echo "12 queued ladders;10 one-thread workers each;maximum6 concurrent ladders (60 cores)."
echo "Each ladder waits for all10 cores and40G memory on one node. Actual concurrency depends on availability."
echo "Task IDs interleave flip/local/uniform by start. This supports a mixed start, not synchronized progress."
echo "50,000 sweeps per replica;burn25,000;no thinning;120h/task;9 sequential adjacent exchange attempts per sweep."
echo "Production is also gated by a matching passed cluster audit, independently of scheduler holds."
qstat -u "$USER"
