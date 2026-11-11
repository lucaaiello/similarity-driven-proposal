#!/bin/bash
# Optional new production run only. Execute from a NEW cluster workspace.
set -euo pipefail
root='Real Data Analysis'
test -f "$root/Data/Prepared/analysis_plan.rds"
sha256sum -c "$root/Cluster/source_files.sha256"
for script in audit ladders summary; do
    bash -n "$root/Cluster/$script.sbatch"
    if command -v h2 >/dev/null 2>&1; then h2 validate "$root/Cluster/$script.sbatch"; fi
done
for path in "$root/Validation/Results/Engine" "$root/Validation/Results/Analysis" "$root/Results/New_Run"; do
    if [ -e "$path" ]; then echo "STOP: existing run/check output: $path. Inspect before resubmitting."; exit 1; fi
done
mkdir -p "$root/Cluster/Logs"
record="$root/Cluster/submitted_jobs.txt"
(set -o noclobber; printf 'THREE METHODS, FOUR FRESH STARTS EACH\n' > "$record") || {
    echo 'Submission record exists. Do not submit duplicate jobs.'; exit 1;
}
audit=$(sbatch --parsable "$root/Cluster/audit.sbatch"); aid=${audit%%;*}
[[ "$aid" =~ ^[0-9]+$ ]]; printf 'audit %s\n' "$aid" >> "$record"
ladders=$(sbatch --parsable --dependency="afterok:$aid" "$root/Cluster/ladders.sbatch"); lid=${ladders%%;*}
[[ "$lid" =~ ^[0-9]+$ ]]; printf 'ladders %s\n' "$lid" >> "$record"
summary=$(sbatch --parsable --dependency="afterok:$lid" "$root/Cluster/summary.sbatch"); sid=${summary%%;*}
[[ "$sid" =~ ^[0-9]+$ ]]; printf 'summary %s\n' "$sid" >> "$record"
cat "$record"
echo '12 ladders; 10 one-thread workers each; up to 6 concurrent. Audit gates and afterok dependencies are required.'
squeue --me
