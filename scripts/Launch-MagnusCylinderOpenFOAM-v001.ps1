param(
 [Parameter(Mandatory=$true)][string]$Authorization,
 [Parameter(Mandatory=$true)][string]$ExpectedPreparationSha256,
 [Parameter(Mandatory=$true)][string]$ExpectedWrapperSha256,
 [Parameter(Mandatory=$true)][string]$ExpectedCommonSha256,
 [Parameter(Mandatory=$true)][string]$ExpectedLauncherSha256,
 [Parameter(Mandatory=$true)][ValidateSet('PrepareOnly','PilotRun')][string]$Mode
)
$ErrorActionPreference='Stop'
. 'E:\CFD-Learning-Lab\local\Use-Project.ps1'
$taskPrivate='E:\CFD-Learning-Lab\local\magnus-cylinder-v001'
$taskCase='E:\CFD-Learning-Lab\native\runs\2026-11-02-magnus-cylinder-v001\openfoam-pilot-v001'
$source='E:\CFD-Learning-Lab\repo\cfd-learning-lab\scripts\prepare_openfoam_magnus_cylinder_v001.py'
$wrapper=Join-Path $taskPrivate 'run-openfoam-v001.sh'
$common='E:\CFD-Learning-Lab\repo\cfd-learning-lab\scripts\magnus_openfoam_common_v001.py'
$authFull=(Resolve-Path -LiteralPath $Authorization).Path
if(!$authFull.StartsWith($taskPrivate+'\',[StringComparison]::OrdinalIgnoreCase)){throw 'Authorization must be in dedicated Magnus private folder'}
$auth=Get-Content -LiteralPath $authFull -Raw|ConvertFrom-Json
if($auth.approved -ne $true -or $auth.approved_by -cne 'root' -or $auth.action -cne 'run_magnus_cylinder_openfoam_v001' -or $auth.allowed_mode -cne $Mode){throw 'Separate exact root phase admission required'}
foreach($pair in @(@($source,$ExpectedPreparationSha256),@($wrapper,$ExpectedWrapperSha256),@($common,$ExpectedCommonSha256),@($PSCommandPath,$ExpectedLauncherSha256))){if((Get-FileHash -LiteralPath $pair[0] -Algorithm SHA256).Hash.ToLowerInvariant() -cne $pair[1]){throw ('Frozen source changed '+$pair[0])}}
if($auth.common_sha256 -cne $ExpectedCommonSha256 -or $auth.preparation_sha256 -cne $ExpectedPreparationSha256 -or $auth.wrapper_sha256 -cne $ExpectedWrapperSha256 -or $auth.launcher_sha256 -cne $ExpectedLauncherSha256){throw 'Authorization must bind all 4 reviewed source hashes'}
if(@(Get-Process fluent*,cortex* -ErrorAction SilentlyContinue).Count){throw 'Fluent must be closed before any WSL phase'}
$ram=[int64](Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory*1024;$disk=[int64](Get-PSDrive E).Free
if($ram -lt 4GB -or $disk -lt 10GB){throw 'Require4GiBavailableRAM/10GiBE before native phase'}
if(!(Test-Path -LiteralPath (Join-Path $taskCase 'provenance\prepared-input-sha256.json'))){throw 'Actual reviewed Windows preparation required'}
if($Mode -eq 'PilotRun'){$receipt=Join-Path $taskCase 'provenance\preflight-validation.json';if(!(Test-Path -LiteralPath $receipt) -or (Get-FileHash -LiteralPath $receipt -Algorithm SHA256).Hash.ToLowerInvariant() -cne $auth.preflight_receipt_sha256){throw 'Exact reviewed actual preflight receipt required'}}
if((Get-FileHash -LiteralPath (Join-Path $taskCase 'provenance\prepared-input-sha256.json') -Algorithm SHA256).Hash.ToLowerInvariant() -cne $auth.prepared_input_manifest_sha256){throw 'Exact reviewed Windows prepared manifest required'}
$nativeMode=if($Mode -eq 'PrepareOnly'){'--prepare-only'}else{'--pilot-run'}
$stamp=Get-Date -Format 'yyyyMMdd-HHmmss-fff';$prefix=Join-Path $taskPrivate ('openfoam-'+$Mode.ToLowerInvariant()+'-v001-'+$stamp)
foreach($suffix in @('-preflight.json','-launch.json','.stdout.txt','.stderr.txt')){if(Test-Path -LiteralPath ($prefix+$suffix)){throw 'Preserve existing launch output'}}
$authLinux='/mnt/e/'+$authFull.Substring(3).Replace('\','/')
@{time=(Get-Date -Format o);mode=$Mode;available_ram_bytes=$ram;free_E_bytes=$disk;mpi_ranks=2;fluent_active=$false;preparation_sha256=$ExpectedPreparationSha256;wrapper_sha256=$ExpectedWrapperSha256;common_sha256=$ExpectedCommonSha256;launcher_sha256=$ExpectedLauncherSha256;authorization_sha256=(Get-FileHash -LiteralPath $authFull -Algorithm SHA256).Hash.ToLowerInvariant();production_admitted=$false}|ConvertTo-Json|Set-Content -LiteralPath ($prefix+'-preflight.json') -Encoding utf8
$child=Start-Process -FilePath 'wsl.exe' -ArgumentList @('-d','CFD-Foundation14','--exec','bash','/mnt/e/CFD-Learning-Lab/local/magnus-cylinder-v001/run-openfoam-v001.sh',$nativeMode,$authLinux) -WorkingDirectory $taskPrivate -WindowStyle Hidden -RedirectStandardOutput ($prefix+'.stdout.txt') -RedirectStandardError ($prefix+'.stderr.txt') -PassThru
$record=@{pid=$child.Id;time=(Get-Date -Format o);mode=$Mode;hidden=$true;mpi_ranks=2;stdout=$prefix+'.stdout.txt';stderr=$prefix+'.stderr.txt';production_admitted=$false}
$record|ConvertTo-Json|Set-Content -LiteralPath ($prefix+'-launch.json') -Encoding utf8
$record|ConvertTo-Json -Compress
