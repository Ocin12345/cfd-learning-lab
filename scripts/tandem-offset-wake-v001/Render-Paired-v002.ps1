param(
    [Parameter(Mandatory=$true)][string]$FluentAudit,
    [Parameter(Mandatory=$true)][string]$FluentManifest,
    [Parameter(Mandatory=$true)][string]$FluentFields,
    [Parameter(Mandatory=$true)][ValidatePattern('^[a-f0-9]{64}$')][string]$FluentFieldsSha256,
    [Parameter(Mandatory=$true)][string]$OpenfoamAudit,
    [Parameter(Mandatory=$true)][string]$OpenfoamManifest,
    [Parameter(Mandatory=$true)][string]$OpenfoamFields,
    [Parameter(Mandatory=$true)][ValidatePattern('^[a-f0-9]{64}$')][string]$OpenfoamFieldsSha256
)
$ErrorActionPreference = 'Stop'
. 'E:/CFD-Learning-Lab/local/Use-Project.ps1'
$plotSource = 'E:/CFD-Learning-Lab/repo/cfd-learning-lab/scripts/plot_tandem_offset_wake_native_v002.py'
$plotDigest = '420820b2991c90f7fcb72ad578a6604620c124edaa20f529f7978c2c13565bce'
if ((Get-FileHash -LiteralPath $plotSource -Algorithm SHA256).Hash.ToLowerInvariant() -ne $plotDigest) { throw 'Reviewed plotter changed' }
$reader = 'E:/CFD-Learning-Lab/repo/cfd-learning-lab/scripts/tandem_offset_wake_fluent_common_v001.py'
if ((Get-FileHash -LiteralPath $reader -Algorithm SHA256).Hash.ToLowerInvariant() -ne '46a5b53aad554a6072a7407a78dc21a9c6a2c687b89b5295cbb85d6b9b3d832a') { throw 'Reviewed native reader changed' }
$exportPath = 'E:/CFD-Learning-Lab/exports/figures/2026-11-02-tandem-offset-wake-paired-v001'
if (Test-Path -LiteralPath $exportPath) { throw 'Fresh paired output directory required' }
$plotArgs = @(
    '-B', $plotSource, '--phase', 'full', '--source-sha256', $plotDigest,
    '--geometry', 'E:/CFD-Learning-Lab/native/runs/2026-11-02-tandem-offset-wake-v001/mesh/geometry.npz',
    '--targets', 'E:/CFD-Learning-Lab/local/tandem-offset-wake-v001/TARGETS.md',
    '--out', $exportPath, '--fluent-reader', $reader,
    '--fluent-audit', $FluentAudit, '--fluent-manifest', $FluentManifest,
    '--fluent-fields', $FluentFields, '--fluent-fields-sha256', $FluentFieldsSha256,
    '--openfoam-audit', $OpenfoamAudit, '--openfoam-manifest', $OpenfoamManifest,
    '--openfoam-fields', $OpenfoamFields, '--openfoam-fields-sha256', $OpenfoamFieldsSha256
)
# Preparation only until root selects actual full audits and admits execution.
# Python performs live >=3 GiB admission, >=2.5 GiB work/10 GiB E: floor checks.
& 'E:/CFD-Learning-Lab/environments/fluent/Scripts/python.exe' @plotArgs
if ($LASTEXITCODE -ne 0) { throw "Actual paired export failed with exit $LASTEXITCODE; preserve source/output" }
