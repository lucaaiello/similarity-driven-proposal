param([string]$Output = 'New_simulations/dm_apt50k_slurm_v1_bundle.tgz')
$ErrorActionPreference = 'Stop'
$repo = Split-Path (Split-Path (Split-Path $PSScriptRoot -Parent) -Parent) -Parent
Push-Location $repo
try {
    if (Test-Path -LiteralPath $Output) { throw "Preserve existing bundle: $Output" }
    $root = 'New_simulations/dm_apt50k_fresh_v1'
    foreach ($line in Get-Content -LiteralPath "$root/inputs/bundle_files.sha256") {
        if ($line -notmatch '^([a-f0-9]{64})  (.+)$') { throw 'Malformed original manifest' }
        $expected = $Matches[1]
        $file = $Matches[2]
        if ((Get-FileHash -LiteralPath $file -Algorithm SHA256).Hash.ToLowerInvariant() -ne $expected) {
            throw "Original source/input changed: $file"
        }
    }
    $extra = @('environment.sh','audit.sbatch','ladders.sbatch','summary.sbatch','submit.sh','README.txt','prepare_bundle.ps1') |
        ForEach-Object { "$root/slurm/$_" }
    $lines = foreach ($file in $extra) {
        $hash = (Get-FileHash -LiteralPath $file -Algorithm SHA256).Hash.ToLowerInvariant()
        "$hash  $file"
    }
    $manifest = "$root/slurm/files.sha256"
    [IO.File]::WriteAllText((Join-Path $repo $manifest), (($lines -join "`n") + "`n"), [Text.Encoding]::ASCII)
    $files = @(Get-Content -LiteralPath 'New_simulations/dm_apt50k_fresh_v1_bundle.tgz.files.txt') + $extra + $manifest
    if (@($files | Select-Object -Unique).Count -ne $files.Count) { throw 'Duplicate archive members' }
    [IO.File]::WriteAllText((Join-Path $repo "$Output.files.txt"), (($files -join "`n") + "`n"), [Text.Encoding]::ASCII)
    & tar.exe -czf $Output -T "$Output.files.txt"
    if ($LASTEXITCODE -ne 0) { throw 'Bundle creation failed' }
    $members = @(& tar.exe -tzf $Output)
    if ($LASTEXITCODE -ne 0 -or @(Compare-Object $files $members).Count -ne 0) { throw 'Archive members differ' }
    $hash = (Get-FileHash -LiteralPath $Output -Algorithm SHA256).Hash.ToLowerInvariant()
    [IO.File]::WriteAllText((Join-Path $repo "$Output.sha256"), "$hash  $([IO.Path]::GetFileName($Output))`n", [Text.Encoding]::ASCII)
    Write-Host "SLURM BUNDLE READY: $($files.Count) explicit files; original source/input manifest verified."
    Write-Host 'No R, compiler, MCMC, SSH or scheduler was executed.'
    Get-FileHash -LiteralPath $Output -Algorithm SHA256
}
finally { Pop-Location }
