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
$plotSource = 'E:/CFD-Learning-Lab/local/magnus-cylinder-v001/plot_magnus_cylinder_native_v002.py'
$plotDigest = '5c40750ce6a502b4634158ce89ddb4e204786890df478b25b6146bfbb50b4b27'
if ((Get-FileHash -LiteralPath $plotSource -Algorithm SHA256).Hash.ToLowerInvariant() -ne $plotDigest) { throw 'Reviewed plotter changed' }
$fluentReader = 'E:/CFD-Learning-Lab/repo/cfd-learning-lab/scripts/magnus_cylinder_fluent_common_v001.py'
$openfoamReader = 'E:/CFD-Learning-Lab/repo/cfd-learning-lab/scripts/magnus_openfoam_common_v001.py'
if ((Get-FileHash -LiteralPath $fluentReader -Algorithm SHA256).Hash.ToLowerInvariant() -ne '1a61e08a2f3624e6ccd36b3c5f244c67d3cbd89eee4865c2a67d030378905820') { throw 'Reviewed Fluent reader changed' }
if ((Get-FileHash -LiteralPath $openfoamReader -Algorithm SHA256).Hash.ToLowerInvariant() -ne '587991b8159edce429e2133170938bbc7c40ece58a0197f1ef63ba0147c3490f') { throw 'Reviewed OpenFOAM reader changed' }
$exportPath = 'E:/CFD-Learning-Lab/exports/figures/2026-11-02-magnus-cylinder-paired-v001'
if (Test-Path -LiteralPath $exportPath) { throw 'Fresh paired output required' }
$plotArgs = @(
    '-B', $plotSource, '--phase', 'full', '--source-sha256', $plotDigest,
    '--geometry', 'E:/CFD-Learning-Lab/native/runs/2026-11-02-magnus-cylinder-v001/mesh/geometry.npz',
    '--targets', 'E:/CFD-Learning-Lab/local/magnus-cylinder-v001/TARGETS.md',
    '--plan', 'E:/CFD-Learning-Lab/local/magnus-cylinder-v001/EXECUTION-PLAN-v001.md',
    '--setup', 'E:/CFD-Learning-Lab/native/runs/2026-11-02-magnus-cylinder-v001/mesh/setup.json',
    '--out', $exportPath, '--fluent-reader', $fluentReader, '--openfoam-reader', $openfoamReader,
    '--fluent-audit', $FluentAudit, '--fluent-manifest', $FluentManifest,
    '--fluent-fields', $FluentFields, '--fluent-fields-sha256', $FluentFieldsSha256,
    '--openfoam-audit', $OpenfoamAudit, '--openfoam-manifest', $OpenfoamManifest,
    '--openfoam-fields', $OpenfoamFields, '--openfoam-fields-sha256', $OpenfoamFieldsSha256
)
# Source preparation only; physical rendering requires root-selected actual full audits.
# Python checks live >=3 GiB RAM admission, >=2.5 GiB during work and >=10 GiB E:.
& 'E:/CFD-Learning-Lab/environments/fluent/Scripts/python.exe' @plotArgs
if ($LASTEXITCODE -ne 0) { throw "Actual paired export failed with exit $LASTEXITCODE; preserve source/output" }
