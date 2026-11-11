param(
    [ValidateSet('parse','compile','fixed','quick','full')][string]$Mode = 'full',
    [ValidateSet(1,2,5,10)][int]$Workers = 2,
    [string]$Rscript = 'C:\Program Files\R\R-4.4.1\bin\Rscript.exe'
)
$ErrorActionPreference = 'Stop'
$repo = Split-Path (Split-Path $PSScriptRoot -Parent) -Parent
Push-Location $repo
try {
    Write-Host 'SEPARATE COMPUTATIONAL COPY. Existing scientific results and sampler sources stay unchanged.'
    if ($Mode -in @('quick','full')) {
        Write-Host 'This verification DOES run bounded short synthetic MCMC;full mode adds six paired real-data single-step fixtures.'
    }
    Write-Host 'No production run or cluster submission is performed.'
    & $Rscript 'New_simulations/dm_apt_runtime_optimization/verify.R' "--mode=$Mode" "--workers=$Workers"
    if ($LASTEXITCODE -ne 0) { throw "Persistent tempering verification failed with exit $LASTEXITCODE. Inspect the new audit output." }
}
finally { Pop-Location }
