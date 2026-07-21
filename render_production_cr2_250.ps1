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
    [string]$CondaExe,
    [string]$RenderEnvironment = 'pymol',
    [string]$FfmpegEnvironment = 'venus_qmmm',
    [string]$FfmpegExe,
    [int]$ChunkSize = 25,
    [int]$Fps = 12
)

$ErrorActionPreference = 'Stop'
$forward = @{
    Trajectory = $Trajectory
    TopologyPdb = $TopologyPdb
    DensityDx = $DensityDx
    RenderRoot = $RenderRoot
    OutputStem = 'tandem_nvt_production_cr2_250_steom_density'
    RenderEnvironment = $RenderEnvironment
    FfmpegEnvironment = $FfmpegEnvironment
    FrameTotal = 250
    FrameDigits = 3
    ChunkSize = $ChunkSize
    Fps = $Fps
}
if (-not [string]::IsNullOrWhiteSpace($MonomerPdb)) { $forward.MonomerPdb = $MonomerPdb }
if (-not [string]::IsNullOrWhiteSpace($OutputVideo)) { $forward.OutputVideo = $OutputVideo }
if (-not [string]::IsNullOrWhiteSpace($CondaExe)) { $forward.CondaExe = $CondaExe }
if (-not [string]::IsNullOrWhiteSpace($FfmpegExe)) { $forward.FfmpegExe = $FfmpegExe }

& (Join-Path $PSScriptRoot 'render_production_cr2_1000.ps1') @forward
