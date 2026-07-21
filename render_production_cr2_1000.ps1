[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$Trajectory,

    [Parameter(Mandatory = $true)]
    [string]$TopologyPdb,

    [Parameter(Mandatory = $true)]
    [string]$DensityDx,

    [Parameter(Mandatory = $true)]
    [string]$RenderRoot,

    [string]$MonomerPdb,
    [string]$OutputVideo,
    [string]$OutputStem = 'tandem_nvt_production_cr2_1000_steom_density',
    [string]$CondaExe,
    [string]$RenderEnvironment = 'pymol',
    [string]$FfmpegEnvironment = 'venus_qmmm',
    [string]$FfmpegExe,
    [int]$FrameTotal = 1000,
    [int]$FrameDigits = 4,
    [int]$ChunkSize = 25,
    [int]$Fps = 12
)

$ErrorActionPreference = 'Stop'
$repoRoot = $PSScriptRoot
if ([string]::IsNullOrWhiteSpace($MonomerPdb)) {
    $MonomerPdb = Join-Path $repoRoot 'tc_simple_old\classical_relaxed.pdb'
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
foreach ($inputPath in @($Trajectory, $TopologyPdb, $DensityDx, $MonomerPdb)) {
    if (-not (Test-Path -LiteralPath $inputPath -PathType Leaf)) {
        throw "Required input does not exist: $inputPath"
    }
}
if ($FrameTotal -le 0 -or $FrameDigits -le 0 -or $ChunkSize -le 0 -or $Fps -le 0) {
    throw 'Frame count, digit count, chunk size, and frame rate must be positive.'
}

$RenderRoot = [System.IO.Path]::GetFullPath($RenderRoot)
$scratchVideo = Join-Path $RenderRoot ($OutputStem + '.mp4')
if ([string]::IsNullOrWhiteSpace($OutputVideo)) {
    $OutputVideo = $scratchVideo
}
else {
    $OutputVideo = [System.IO.Path]::GetFullPath($OutputVideo)
}
$frames = Join-Path $RenderRoot ($OutputStem + '_frames')
$logs = Join-Path $RenderRoot 'logs'
$outputParent = Split-Path -Parent $OutputVideo
New-Item -ItemType Directory -Path $RenderRoot, $logs, $outputParent -Force | Out-Null

$renderer = Join-Path $repoRoot 'render_nvt_movie.py'
$frameName = 'frm_{0:D' + $FrameDigits + '}.png'
for ($start = 1; $start -le $FrameTotal; $start += $ChunkSize) {
    $end = [Math]::Min($FrameTotal, $start + $ChunkSize - 1)
    $expected = $start..$end | ForEach-Object { Join-Path $frames ($frameName -f $_) }
    $complete = ($expected | Where-Object { -not (Test-Path -LiteralPath $_) }).Count -eq 0
    if ($complete) {
        Write-Output ("[render] chunk $start-$end already complete")
        continue
    }

    $log = Join-Path $logs ("chunk_{0}_{1}.log" -f $start, $end)
    $arguments = @(
        'run', '--no-capture-output', '-n', $RenderEnvironment, 'pymol', '-cq', $renderer, '--',
        '--traj', (Resolve-Path -LiteralPath $Trajectory).Path,
        '--topology-pdb', (Resolve-Path -LiteralPath $TopologyPdb).Path,
        '--out', $scratchVideo,
        '--rep', 'density',
        '--density-dx', (Resolve-Path -LiteralPath $DensityDx).Path,
        '--monomer-pdb', (Resolve-Path -LiteralPath $MonomerPdb).Path,
        '--width', '1920', '--height', '1440', '--fps', [string]$Fps,
        '--stride', '1', '--start-frame', [string]$start, '--end-frame', [string]$end,
        '--iso-level', '0.010', '--iso-transparency', '0.18',
        '--rotate-deg', '90', '--keep-frames', '--frames-only'
    )
    & $CondaExe @arguments *> $log
    if ($LASTEXITCODE -ne 0) {
        throw "PyMOL render failed for $start-$end; see $log"
    }
    $missing = $expected | Where-Object { -not (Test-Path -LiteralPath $_) }
    if ($missing) {
        throw "PyMOL returned success but chunk $start-$end is incomplete"
    }
    Write-Output ("[render] completed chunk $start-$end")
}

$frameCount = (Get-ChildItem -LiteralPath $frames -Filter 'frm_*.png').Count
if ($frameCount -ne $FrameTotal) {
    throw "Expected $FrameTotal PNG frames, found $frameCount"
}

$ffmpegArguments = @(
    '-y', '-framerate', [string]$Fps, '-start_number', '1',
    '-i', (Join-Path $frames ("frm_%0${FrameDigits}d.png")),
    '-c:v', 'libx264', '-preset', 'slow', '-crf', '18', '-pix_fmt', 'yuv420p',
    '-movflags', '+faststart', $scratchVideo
)
if ([string]::IsNullOrWhiteSpace($FfmpegExe)) {
    $ffmpegCommand = Get-Command ffmpeg -ErrorAction SilentlyContinue
    if ($null -ne $ffmpegCommand) {
        $FfmpegExe = $ffmpegCommand.Source
    }
}
if (-not [string]::IsNullOrWhiteSpace($FfmpegExe)) {
    & $FfmpegExe @ffmpegArguments
}
else {
    & $CondaExe run --no-capture-output -n $FfmpegEnvironment ffmpeg @ffmpegArguments
}
if ($LASTEXITCODE -ne 0) {
    throw 'ffmpeg encoding failed'
}

if ([System.IO.Path]::GetFullPath($scratchVideo) -ne [System.IO.Path]::GetFullPath($OutputVideo)) {
    Copy-Item -LiteralPath $scratchVideo -Destination $OutputVideo -Force
}
Write-Output "[render] encoded $frameCount frames -> $OutputVideo"
