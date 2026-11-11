[CmdletBinding()]
param(
    [string]$Rscript = 'C:\Program Files\R\R-4.4.1\bin\Rscript.exe',
    [string]$OutputRoot = '',
    [switch]$PilotOnly
)
$ErrorActionPreference = 'Stop'
$repoRoot = Split-Path -Parent $PSScriptRoot
$originalDirectory = Get-Location
$originalThreads = @{}
foreach ($threadName in @('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')) {
    $originalThreads[$threadName] = [Environment]::GetEnvironmentVariable($threadName,'Process')
}
$timings = [System.Collections.Generic.List[object]]::new()
function Invoke-RevisionStep {
    param([string]$Name, [string]$Script, [string[]]$RArguments)
    Write-Host ("`n[{0}] {1}" -f (Get-Date -Format 'HH:mm:ss'), $Name)
    $stopwatch = [Diagnostics.Stopwatch]::StartNew()
    # Windows PowerShell otherwise treats harmless native stderr warnings as
    # terminating PowerShell errors. Native exit status is checked explicitly.
    $previousPreference = $ErrorActionPreference
    try {
        $ErrorActionPreference = 'Continue'
        & $Rscript $Script @RArguments 2>&1 | ForEach-Object { $_.ToString() } | Tee-Object -FilePath (Join-Path $OutputRoot ($Name + '.log'))
        $nativeExitCode = $LASTEXITCODE
    } finally {
        $ErrorActionPreference = $previousPreference
        $stopwatch.Stop()
    }
    $timings.Add([pscustomobject]@{Step=$Name;Seconds=$stopwatch.Elapsed.TotalSeconds;ExitCode=$nativeExitCode})
    $timings | Export-Csv -LiteralPath (Join-Path $OutputRoot 'pipeline_timings.csv') -NoTypeInformation
    if ($nativeExitCode -ne 0) { throw "$Name failed (exit $nativeExitCode). No subsequent experiment was started." }
}
try {
    Set-Location -LiteralPath $repoRoot
    if (-not (Test-Path -LiteralPath $Rscript -PathType Leaf)) { throw "Rscript not found: $Rscript" }
    if (-not $OutputRoot) {
        $OutputRoot = Join-Path $repoRoot ('New_simulations\outputs\matched_gaussian_' + (Get-Date -Format 'yyyyMMdd_HHmmss'))
    } else {
        $OutputRoot = [IO.Path]::GetFullPath($OutputRoot)
    }
    if (Test-Path -LiteralPath $OutputRoot) { throw "Choose a NEW output folder; refusing to overwrite: $OutputRoot" }
    New-Item -ItemType Directory -Path $OutputRoot | Out-Null
    foreach ($threadName in $originalThreads.Keys) {
        [Environment]::SetEnvironmentVariable($threadName,'1','Process')
    }
    Write-Host "Results: $OutputRoot"
    $kernelAudit = Join-Path $OutputRoot 'audit_kernel'
    $rankAudit = Join-Path $OutputRoot 'audit_rank_one'
    $swapAudit = Join-Path $OutputRoot 'audit_swap'
    $marginalAudit = Join-Path $OutputRoot 'audit_marginal'
    Invoke-RevisionStep '01_kernel_audit' 'New_simulations/run_matched_gaussian_revision_audit.R' @("--output=$kernelAudit")
    Invoke-RevisionStep '02_rank_one_audit' 'New_simulations/run_production_rank_one_equivalence_audit.R' @("--output=$rankAudit")
    Invoke-RevisionStep '03_swap_audit' 'New_simulations/run_local_move_rank_one_equivalence_audit.R' @("--output=$swapAudit")
    Invoke-RevisionStep '04_marginal_audit' 'New_simulations/run_marginal_score_equivalence_audit.R' @("--output=$marginalAudit")
    $gate = "--revision-audit=$kernelAudit"
    $rankGate = "--production-audit=$(Join-Path $rankAudit 'summary.csv')"
    $marginalGate = "--marginal-audit=$(Join-Path $marginalAudit 'summary.csv')"
    $swapGate = "--audit=$(Join-Path $swapAudit 'summary.csv')"
    Invoke-RevisionStep '05_exact_pilot' 'New_simulations/run_production_adaptive_exact_benchmark.R' @('--mode=quick','--design=grid','--method-set=all','--proposal-engine=rank_one','--datasets=1',$gate,"--output=$(Join-Path $OutputRoot 'pilot_exact')")
    Invoke-RevisionStep '06_scaling_pilot' 'New_simulations/run_highdim_all_methods_scaling.R' @('--mode=quick','--dimensions=50,200,1000','--datasets=1','--chains=1','--budget=0.25',$gate,"--audit-dir=$rankAudit","--output=$(Join-Path $OutputRoot 'pilot_scaling')")
    Invoke-RevisionStep '07_swap_pilot' 'New_simulations/run_correlated_surrogate_local_move_ablation.R' @('--mode=quick','--datasets=1',$gate,$swapGate,"--output=$(Join-Path $OutputRoot 'pilot_swap')")
    Invoke-RevisionStep '08_exact_score_pilot' 'New_simulations/run_auxiliary_score_exact_sensitivity.R' @('--mode=quick',$gate,$rankGate,$marginalGate,"--output=$(Join-Path $OutputRoot 'pilot_exact_score')")
    Invoke-RevisionStep '09_correlated_score_pilot' 'New_simulations/run_auxiliary_score_correlated_robustness.R' @('--mode=quick',$gate,$rankGate,$marginalGate,"--output=$(Join-Path $OutputRoot 'pilot_correlated_score')")
    if (-not $PilotOnly) {
        Invoke-RevisionStep '10_exact_full' 'New_simulations/run_production_adaptive_exact_benchmark.R' @('--mode=full','--design=grid','--method-set=all','--proposal-engine=rank_one','--datasets=20',$gate,"--output=$(Join-Path $OutputRoot 'exact')")
        Invoke-RevisionStep '11_scaling_full' 'New_simulations/run_highdim_all_methods_scaling.R' @('--mode=full',$gate,"--audit-dir=$rankAudit","--output=$(Join-Path $OutputRoot 'scaling')")
        Invoke-RevisionStep '12_swap_full' 'New_simulations/run_correlated_surrogate_local_move_ablation.R' @('--mode=full',$gate,$swapGate,"--output=$(Join-Path $OutputRoot 'swap')")
        Invoke-RevisionStep '13_exact_score_full' 'New_simulations/run_auxiliary_score_exact_sensitivity.R' @('--mode=full',$gate,$rankGate,$marginalGate,"--output=$(Join-Path $OutputRoot 'exact_score')")
        Invoke-RevisionStep '14_correlated_score_full' 'New_simulations/run_auxiliary_score_correlated_robustness.R' @('--mode=full',$gate,$rankGate,$marginalGate,"--output=$(Join-Path $OutputRoot 'correlated_score')")
        Invoke-RevisionStep '15_figures_and_tables' 'New_simulations/make_manuscript_simulation_figures.R' @('--figure=gaussian',"--gaussian-root=$OutputRoot")
    }
    Write-Host "`nCompleted. Results and per-stage logs: $OutputRoot"
} finally {
    foreach ($threadName in $originalThreads.Keys) {
        [Environment]::SetEnvironmentVariable($threadName,$originalThreads[$threadName],'Process')
    }
    Set-Location -LiteralPath $originalDirectory
}
