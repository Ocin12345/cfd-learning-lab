param(
 [Parameter(Mandatory=$true)][string]$Authorization,
 [Parameter(Mandatory=$true)][string]$ExpectedPreparationSha256,
 [Parameter(Mandatory=$true)][string]$ExpectedWrapperSha256,
 [Parameter(Mandatory=$true)][string]$ExpectedLauncherSha256,
 [Parameter(Mandatory=$true)][ValidateSet('ProductionRun')][string]$Mode
)
$ErrorActionPreference='Stop'
. 'E:\CFD-Learning-Lab\local\Use-Project.ps1'
$taskPrivate='E:\CFD-Learning-Lab\local\tandem-offset-wake-v001'
$taskCase='E:\CFD-Learning-Lab\native\runs\2026-11-02-tandem-offset-wake-v001\openfoam-production-v001'
$source='E:\CFD-Learning-Lab\repo\cfd-learning-lab\scripts\prepare_openfoam_tandem_offset_wake_production_v001.py'
$wrapper=Join-Path $taskPrivate 'run-openfoam-production-v001.sh'
$authFull=(Resolve-Path -LiteralPath $Authorization).Path
if(!$authFull.StartsWith($taskPrivate+'\',[StringComparison]::OrdinalIgnoreCase)){throw 'Authorization must be in dedicated tandem private folder'}
$auth=Get-Content -LiteralPath $authFull -Raw|ConvertFrom-Json
if($auth.approved -ne $true -or $auth.approved_by -cne 'root' -or $auth.action -cne 'run_tandem_offset_wake_production_openfoam_v001' -or $auth.allowed_mode -cne $Mode){throw 'Separate exact root phase admission required'}
foreach($pair in @(@($source,$ExpectedPreparationSha256),@($wrapper,$ExpectedWrapperSha256),@($PSCommandPath,$ExpectedLauncherSha256))){if((Get-FileHash -LiteralPath $pair[0] -Algorithm SHA256).Hash.ToLowerInvariant() -cne $pair[1]){throw ('Frozen source changed '+$pair[0])}}
if($auth.preparation_sha256 -cne $ExpectedPreparationSha256 -or $auth.wrapper_sha256 -cne $ExpectedWrapperSha256 -or $auth.launcher_sha256 -cne $ExpectedLauncherSha256){throw 'Authorization must bind all 3 reviewed production source hashes'}
if(@(Get-Process fluent*,cortex* -ErrorAction SilentlyContinue).Count){throw 'Fluent must be closed before any WSL phase'}
$ram=[int64](Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory*1024;$disk=[int64](Get-PSDrive E).Free
if($ram -lt 4GB -or $disk -lt 10GB){throw 'Require4GiBavailableRAM/10GiBE before native phase'}
if(!(Test-Path -LiteralPath (Join-Path $taskCase 'provenance\prepared-primary-sha256.json'))){throw 'Actual reviewed Windows preparation required'}

foreach($pair in @(@('prepared_primary_manifest_sha256','provenance\prepared-primary-sha256.json'),@('inherited_accepted_pilot_receipt_sha256','provenance\inherited-accepted-pilot.json'),@('accepted_pilot_selection_sha256','provenance\accepted-pilot-selection.json'))){if((Get-FileHash -LiteralPath (Join-Path $taskCase $pair[1]) -Algorithm SHA256).Hash.ToLowerInvariant() -cne $auth.($pair[0])){throw ('Exact reviewed accepted pilot copy required '+$pair[1])}}
$nativeMode='--production-run'
$stamp=Get-Date -Format 'yyyyMMdd-HHmmss-fff';$prefix=Join-Path $taskPrivate ('openfoam-'+$Mode.ToLowerInvariant()+'-v001-'+$stamp)
foreach($suffix in @('-preflight.json','-launch.json','.stdout.txt','.stderr.txt')){if(Test-Path -LiteralPath ($prefix+$suffix)){throw 'Preserve existing launch output'}}
$authLinux='/mnt/e/'+$authFull.Substring(3).Replace('\','/')
@{time=(Get-Date -Format o);mode=$Mode;available_ram_bytes=$ram;free_E_bytes=$disk;mpi_ranks=2;fluent_active=$false;preparation_sha256=$ExpectedPreparationSha256;wrapper_sha256=$ExpectedWrapperSha256;launcher_sha256=$ExpectedLauncherSha256;authorization_sha256=(Get-FileHash -LiteralPath $authFull -Algorithm SHA256).Hash.ToLowerInvariant();production_admitted=$true}|ConvertTo-Json|Set-Content -LiteralPath ($prefix+'-preflight.json') -Encoding utf8
$child=Start-Process -FilePath 'wsl.exe' -ArgumentList @('-d','CFD-Foundation14','--exec','bash','/mnt/e/CFD-Learning-Lab/local/tandem-offset-wake-v001/run-openfoam-production-v001.sh',$nativeMode,$authLinux) -WorkingDirectory $taskPrivate -WindowStyle Hidden -RedirectStandardOutput ($prefix+'.stdout.txt') -RedirectStandardError ($prefix+'.stderr.txt') -PassThru
$record=@{pid=$child.Id;time=(Get-Date -Format o);mode=$Mode;hidden=$true;mpi_ranks=2;stdout=$prefix+'.stdout.txt';stderr=$prefix+'.stderr.txt';production_admitted=$true}
$record|ConvertTo-Json|Set-Content -LiteralPath ($prefix+'-launch.json') -Encoding utf8
$record|ConvertTo-Json -Compress
