param(
 [Parameter(Mandatory=$true)][string]$Authorization,
 [Parameter(Mandatory=$true)][string]$Out,
 [Parameter(Mandatory=$true)][ValidatePattern('^[0-9a-f]{64}$')][string]$SourceFreezeSha256,
 [switch]$PreflightOnly
)
$ErrorActionPreference='Stop'
. 'E:\CFD-Learning-Lab\local\Use-Project.ps1'
$pilotPrivate='E:\CFD-Learning-Lab\local\tandem-offset-wake-v001'
$pilotFreeze=Join-Path $pilotPrivate 'fluent-production-source-freeze-v001.json'
$pilotRunner='E:\CFD-Learning-Lab\repo\cfd-learning-lab\scripts\tandem_offset_wake_fluent_production_v001.py'
$pilotPython='E:\CFD-Learning-Lab\environments\fluent\Scripts\python.exe'
function Assert-FrozenSources {
 if((Get-FileHash -LiteralPath $pilotFreeze -Algorithm SHA256).Hash.ToLowerInvariant() -ne $SourceFreezeSha256){throw 'Root-pinned source freeze differs'}
 $pilotData=Get-Content -LiteralPath $pilotFreeze -Raw | ConvertFrom-Json
 foreach($pilotRow in $pilotData.files){
  if((Get-FileHash -LiteralPath $pilotRow.path -Algorithm SHA256).Hash.ToLowerInvariant() -ne $pilotRow.sha256){throw ('Frozen source differs: '+$pilotRow.path)}
 }
}
function Assert-NoOtherCFD {
 $pilotCFD=@(Get-CimInstance Win32_Process | Where-Object {
  $_.Name -match '^(fluent|cortex|fl\d+|foamRun|simpleFoam|pimpleFoam)(\.exe)?$' -or
  ($_.Name -match '^wsl(\.exe)?$' -and $_.CommandLine -match '(foamRun|simpleFoam|pimpleFoam)')
 })
 if($pilotCFD.Count -ne 0){throw 'Another actual CFD solver is active; sequential launches required'}
}
Assert-FrozenSources
Assert-NoOtherCFD
$pilotAuth=(Resolve-Path -LiteralPath $Authorization).Path
$pilotOut=[System.IO.Path]::GetFullPath($Out)
$pilotPreflight=& $pilotPython -B $pilotRunner --authorization $pilotAuth --out $pilotOut --prepare-check
if($LASTEXITCODE -ne 0){throw 'Pilot root/source/resource read-only preflight refused'}
if($PreflightOnly){$pilotPreflight;return}
Assert-FrozenSources
Assert-NoOtherCFD
$pilotStamp=Get-Date -Format 'yyyyMMdd-HHmmss-fff'
$pilotStdout=Join-Path $pilotPrivate ('fluent-production-'+$pilotStamp+'.stdout.txt')
$pilotStderr=Join-Path $pilotPrivate ('fluent-production-'+$pilotStamp+'.stderr.txt')
$pilotArgs=@('-B',('"'+$pilotRunner+'"'),'--authorization',('"'+$pilotAuth+'"'),'--out',('"'+$pilotOut+'"'))
$pilotProcess=Start-Process -FilePath $pilotPython -ArgumentList $pilotArgs -WorkingDirectory $pilotPrivate -WindowStyle Hidden -PassThru -RedirectStandardOutput $pilotStdout -RedirectStandardError $pilotStderr
$pilotReceipt=[ordered]@{action='tandem_offset_wake_fluent_production_v001';phase='production-795-additional-steps';python_pid=$pilotProcess.Id;out=$pilotOut;authorization=$pilotAuth;authorization_sha256=(Get-FileHash -LiteralPath $pilotAuth).Hash.ToLowerInvariant();source_freeze_sha256=$SourceFreezeSha256;stdout=$pilotStdout;stderr=$pilotStderr;solver_processes=2;maximum_additional_physical_steps=795;maximum_total_steps=800;maximum_end_time_s=40.;production_authorized=$true;publication_performed=$false}
$pilotReceipt | ConvertTo-Json -Depth 10 | Set-Content -LiteralPath (Join-Path $pilotPrivate ('fluent-production-'+$pilotStamp+'-launch.json')) -Encoding utf8
$pilotReceipt | ConvertTo-Json -Depth 10
