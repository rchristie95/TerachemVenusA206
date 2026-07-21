[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$Transforms,

    [Parameter(Mandatory = $true)]
    [string]$Density,

    [string]$Output,
    [string]$CondaExe,
    [string]$CondaEnvironment = 'venus_qmmm',
    [double]$Epsilon = 1.77,
    [ValidateSet('f32', 'f64')]
    [string]$Precision = 'f32',
    [double]$BinWidth = 0.0
)

$ErrorActionPreference = 'Stop'
$repoRoot = $PSScriptRoot
if ([string]::IsNullOrWhiteSpace($Output)) {
    $Output = Join-Path $repoRoot 'coupling_nvt_production_cr2_1000'
}
if ([string]::IsNullOrWhiteSpace($CondaExe)) {
    if (-not [string]::IsNullOrWhiteSpace($env:CONDA_EXE)) {
        $CondaExe = $env:CONDA_EXE
    }
    else {
        $condaCommand = Get-Command conda -ErrorAction SilentlyContinue
        if ($null -ne $condaCommand) {
            $CondaExe = $condaCommand.Source
        }
        else {
            $commonConda = @(
                (Join-Path $env:USERPROFILE 'miniforge3\Scripts\conda.exe'),
                (Join-Path $env:USERPROFILE 'anaconda3\Scripts\conda.exe'),
                (Join-Path $env:USERPROFILE 'miniconda3\Scripts\conda.exe')
            ) | Where-Object { Test-Path -LiteralPath $_ -PathType Leaf } | Select-Object -First 1
            if ($null -eq $commonConda) {
                throw 'Conda was not found. Activate Conda or pass -CondaExe.'
            }
            $CondaExe = $commonConda
        }
    }
}
foreach ($inputPath in @($Transforms, $Density)) {
    if (-not (Test-Path -LiteralPath $inputPath -PathType Leaf)) {
        throw "Required input does not exist: $inputPath"
    }
}
if ($Epsilon -le 0.0 -or $BinWidth -lt 0.0) {
    throw 'Epsilon must be positive and bin width must be non-negative.'
}

$arguments = @(
    'run', '--no-capture-output', '-n', $CondaEnvironment, 'python',
    (Join-Path $repoRoot 'coupling_dcd_steom.py'),
    '--density', (Resolve-Path -LiteralPath $Density).Path,
    '--transforms', (Resolve-Path -LiteralPath $Transforms).Path,
    '--out', [System.IO.Path]::GetFullPath($Output),
    '--epsilon', [string]$Epsilon,
    '--bin-width', [string]$BinWidth,
    '--precision', $Precision
)

& $CondaExe @arguments
if ($LASTEXITCODE -ne 0) {
    throw "Full-grid coupling calculation failed with exit code $LASTEXITCODE"
}
