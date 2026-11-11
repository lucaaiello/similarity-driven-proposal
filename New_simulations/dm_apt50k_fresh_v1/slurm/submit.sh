#!/bin/bash
# Run only in the NEW Slurm workspace, after holding the legacy jobs.
set -euo pipefail
root=New_simulations/dm_apt50k_fresh_v1
test -f "$root/inputs/plan.rds"
test -f "$root/slurm/LEGACY_HELD.txt" || {
    echo "Hold legacy jobs 15351594 and15351595, then record confirmation in slurm/LEGACY_HELD.txt."; exit 1;
}
sha256sum -c "$root/inputs/bundle_files.sha256"
sha256sum -c "$root/slurm/files.sha256"
mkdir -p "$root/slurm/logs"
for script in audit ladders summary; do
    bash -n "$root/slurm/$script.sbatch"
    h2 validate "$root/slurm/$script.sbatch"
done
for path in "$root/engine_audit" "$root/audit" New_simulations/outputs/dm_apt50k_fresh_v1; do
    if [ -e "$path" ]; then
        echo "STOP: existing audit/production output: $path. Inspect rather than resubmit."; exit 1;
    fi
done
record="$root/slurm/submitted_jobs.txt"
( set -o noclobber; printf 'FRESH50K SLURM: THREE METHODS xFOUR FRESH STARTS\n' > "$record" ) || {
    echo "Submission record already exists. Do not submit duplicates."; exit 1;
}
audit=$(sbatch --parsable "$root/slurm/audit.sbatch"); aid=${audit%%;*}
[[ "$aid" =~ ^[0-9]+$ ]]
printf 'audit %s\n' "$aid" >> "$record"
ladders=$(sbatch --parsable --dependency="afterok:$aid" "$root/slurm/ladders.sbatch"); lid=${ladders%%;*}
[[ "$lid" =~ ^[0-9]+$ ]]
printf 'ladders %s\n' "$lid" >> "$record"
summary=$(sbatch --parsable --dependency="afterok:$lid" "$root/slurm/summary.sbatch"); sid=${summary%%;*}
[[ "$sid" =~ ^[0-9]+$ ]]
printf 'summary %s\n' "$sid" >> "$record"
cat "$record"
echo "12 ladders;10 one-thread workers each;up to6 concurrent;40G and120h per ladder."
echo "Audit success is required by afterok AND the unchanged matching-runtime R audit gate."
echo "If any job fails, dependent jobs may remain pending. Inspect; do not rerun this script."
squeue --me
