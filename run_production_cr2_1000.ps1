[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$WorkDir,

    [string]$Pdb,
    [string]$AmberCr2Prmtop,
    [string]$CondaExe,
    [string]$CondaEnvironment = 'venus_qmmm',
    [ValidateSet('CUDA', 'OpenCL', 'CPU')]
    [string]$Platform = 'OpenCL',
    [int]$NvtSteps = 500000,
    [double]$TimestepFs = 2.0,
    [int]$TrajectoryInterval = 500,
    [int]$CheckpointInterval = 10000,
    [string]$TrajectoryFile = 'tandem_nvt_production_cr2_1000.dcd'
)

$ErrorActionPreference = 'Stop'
$repoRoot = $PSScriptRoot
if ([string]::IsNullOrWhiteSpace($Pdb)) {
    $Pdb = Join-Path $repoRoot 'tandem_dimer_production_cr2.pdb'
}
if ([string]::IsNullOrWhiteSpace($AmberCr2Prmtop)) {
    $AmberCr2Prmtop = Join-Path $repoRoot 'anionic_build\monomer_solv.prmtop'
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

foreach ($inputPath in @($Pdb, $AmberCr2Prmtop)) {
    if (-not (Test-Path -LiteralPath $inputPath -PathType Leaf)) {
        throw "Required input does not exist: $inputPath"
    }
}
if ($NvtSteps -le 0 -or $TimestepFs -le 0 -or $TrajectoryInterval -le 0 -or $CheckpointInterval -le 0) {
    throw 'Step counts, timestep, trajectory interval, and checkpoint interval must be positive.'
}

$runner = Join-Path $repoRoot 'run_nvt.py'
$pythonArgs = @(
    'run', '--no-capture-output', '-n', $CondaEnvironment, 'python', '-u',
    $runner,
    '--pdb', (Resolve-Path -LiteralPath $Pdb).Path,
    '--amber-cr2-prmtop', (Resolve-Path -LiteralPath $AmberCr2Prmtop).Path,
    '--workdir', [System.IO.Path]::GetFullPath($WorkDir),
    '--nvt-steps', [string]$NvtSteps,
    '--timestep-fs', [string]$TimestepFs,
    '--openmm-trajectory-file', $TrajectoryFile,
    '--openmm-trajectory-interval', [string]$TrajectoryInterval,
    '--trajectory-protein-only',
    '--checkpoint-interval', [string]$CheckpointInterval,
    '--platform', $Platform,
    '--no-video'
)

& $CondaExe @pythonArgs
if ($LASTEXITCODE -ne 0) {
    throw "Production NVT failed with exit code $LASTEXITCODE"
}
