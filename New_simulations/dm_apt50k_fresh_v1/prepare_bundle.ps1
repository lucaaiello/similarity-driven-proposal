param(
    [string]$Output = 'New_simulations/dm_apt50k_fresh_v1_bundle.tgz',
    [string]$Rscript = 'C:\Program Files\R\R-4.4.1\bin\Rscript.exe'
)
$ErrorActionPreference = 'Stop'
$repo = Split-Path (Split-Path $PSScriptRoot -Parent) -Parent
Push-Location $repo
try {
    if (Test-Path -LiteralPath $Output) { throw "Preserve existing bundle: $Output. Supply another -Output filename." }
    & $Rscript 'New_simulations/dm_apt50k_fresh_v1/prepare.R'
    if ($LASTEXITCODE -ne 0) { throw 'Fresh input/plan preparation failed' }
    $original = 'New_simulations/dm_adaptive_tempering'
    $copy = 'New_simulations/dm_apt_runtime_optimization'
    $fresh = 'New_simulations/dm_apt50k_fresh_v1'
    $files = @('DM_tempered.cpp','joint_kernel.h','dm_stable_numerics.h',
        'baseline/DM_joint_block.cpp','baseline/joint_kernel.h','baseline/dm_stable_numerics.h',
        'support.R','math.R','reference.R','prepare.R','audit.R','runner.R','run_task.R','summarize.R',
        'inputs/plan.rds','inputs/origins.rds') | ForEach-Object { "$original/$_" }
    $files += @('DM_persistent.cpp','joint_kernel.h','dm_stable_numerics.h',
        'baseline/DM_tempered.cpp','baseline/joint_kernel.h','baseline/dm_stable_numerics.h',
        'baseline/support.R','baseline/math.R','baseline/runner.R','baseline/reference.R',
        'audit_cases.R','support.R','verify.R','native_generation.rds') | ForEach-Object { "$copy/$_" }
    $files += @('support.R','prepare.R','runner.R','audit.R','run_task.R','summarize.R',
        'prepare_bundle.ps1','README.txt','inputs/plan.rds','inputs/fresh.rds',
        'hoffman2/audit.qsub','hoffman2/ladders.qsub','hoffman2/summary.qsub','hoffman2/submit.sh') |
        ForEach-Object { "$fresh/$_" }
    foreach ($file in $files) {
        if (-not (Test-Path -LiteralPath $file -PathType Leaf)) { throw "Missing file: $file" }
    }
    $checksums = foreach ($file in $files) {
        $hash = (Get-FileHash -LiteralPath $file -Algorithm SHA256).Hash.ToLowerInvariant()
        "$hash  $file"
    }
    $manifest = "$fresh/inputs/bundle_files.sha256"
    [IO.File]::WriteAllText((Join-Path $repo $manifest), (($checksums -join "`n") + "`n"), [Text.Encoding]::ASCII)
    $files += $manifest
    [IO.File]::WriteAllText((Join-Path $repo "$Output.files.txt"), (($files -join "`n") + "`n"), [Text.Encoding]::ASCII)
    & tar.exe -czf $Output -T "$Output.files.txt"
    if ($LASTEXITCODE -ne 0) { throw 'Bundle creation failed' }
    $members = @(& tar.exe -tzf $Output)
    if ($LASTEXITCODE -ne 0 -or @(Compare-Object $files $members).Count -ne 0) { throw 'Archive file list differs' }
    $hash = (Get-FileHash -LiteralPath $Output -Algorithm SHA256).Hash.ToLowerInvariant()
    [IO.File]::WriteAllText((Join-Path $repo "$Output.sha256"), "$hash  $([IO.Path]::GetFileName($Output))`n", [Text.Encoding]::ASCII)
    Write-Host "BUNDLE READY: $($files.Count) explicit files. No compilation, optimizer, MCMC or submission."
    Get-FileHash -LiteralPath $Output -Algorithm SHA256
}
finally { Pop-Location }
