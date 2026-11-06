# SOURCE DRAFT. Hidden one-phase OpenFOAM qualification; no production mode.
param(
 [Parameter(Mandatory=$true)][ValidateSet('Control','Pulse')][string]$Phase,
 [Parameter(Mandatory=$true)][string]$Authorization,
 [Parameter(Mandatory=$true)][string]$ExpectedPreparationSha256,
 [Parameter(Mandatory=$true)][string]$ExpectedWrapperSha256,
 [Parameter(Mandatory=$true)][string]$ExpectedLauncherSha256
)
$ErrorActionPreference='Stop'
. 'E:\CFD-Learning-Lab\local\Use-Project.ps1'
$private='E:\CFD-Learning-Lab\local\prospective-tracer-v001'
$prep=Join-Path $private 'prepare_openfoam_tracer_qualification_v004.py'
$wrapper=Join-Path $private 'run-openfoam-tracer-qualification-v004.sh'
$python='E:\CFD-Learning-Lab\environments\fluent\Scripts\python.exe'
$authFull=(Resolve-Path -LiteralPath $Authorization).Path
if(!$authFull.StartsWith($private+'\',[StringComparison]::OrdinalIgnoreCase)){throw 'Authorization must remain in exact candidate private folder'}
$auth=Get-Content -LiteralPath $authFull -Raw | ConvertFrom-Json
if($auth.approved -ne $true -or $auth.approved_by -cne 'root' -or $auth.action -cne 'run_tracer_openfoam_qualification_v004' -or $auth.allowed_phase -cne $Phase -or $auth.actual_steps_authorized -ne 5 -or $auth.final_time_s -ne 0.125 -or $auth.production_admitted -ne $false -or $auth.fluent_sessions_closed -ne $true){throw 'Separate exact five-step root phase authorization required'}
foreach($pair in @(@($prep,$ExpectedPreparationSha256),@($wrapper,$ExpectedWrapperSha256),@($PSCommandPath,$ExpectedLauncherSha256))){
 if($pair[1] -cnotmatch '^[a-f0-9]{64}$' -or (Get-FileHash -LiteralPath $pair[0] -Algorithm SHA256).Hash.ToLowerInvariant() -cne $pair[1]){throw ('Reviewed source changed: '+$pair[0])}
}
if($auth.preparation_sha256 -cne $ExpectedPreparationSha256 -or $auth.wrapper_sha256 -cne $ExpectedWrapperSha256 -or $auth.launcher_sha256 -cne $ExpectedLauncherSha256){throw 'Authorization must bind all three exact reviewed source hashes'}
if(@(Get-Process fluent*,cortex*,fl_mpi*,fl2d*,fl3d* -ErrorAction SilentlyContinue).Count){throw 'Fluent must be closed before any WSL/native preflight'}
$ram=[int64](Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory*1024
$disk=[int64](Get-PSDrive E).Free
if($ram -lt 4GB -or $disk -lt 10GB){throw ('Resource admission held: '+[math]::Round($ram/1GB,3)+' GiB RAM; '+[math]::Round($disk/1GB,3)+' GiB E:')}
$caseName=if($Phase -ceq 'Control'){'openfoam-control-v003'}else{'openfoam-pulse-pilot-v003'}
$casePath=Join-Path 'E:\CFD-Learning-Lab\native\runs\2026-11-02-tracer-tjunction-v001' $caseName
if(Test-Path -LiteralPath $casePath){throw 'Fresh versioned case required; preserve existing attempt'}
& $python -B $prep --phase $Phase --authorization $authFull --verify-sources
if($LASTEXITCODE -ne 0){throw 'Qualified ancestry/source/auth preflight failed before WSL start'}
$stamp=Get-Date -Format 'yyyyMMdd-HHmmss-fff'
$prefix=Join-Path $private ('openfoam-'+$Phase.ToLowerInvariant()+'-v004-'+$stamp)
$authLinux='/mnt/e/'+$authFull.Substring(3).Replace('\','/')
@{time=(Get-Date -Format o);phase=$Phase;available_ram_bytes=$ram;ram_floor_bytes=4GB;free_E_bytes=$disk;disk_floor_bytes=10GB;processes=2;fluent_active=$false;native_flow_qualification_sha256=$auth.native_flow_qualification_sha256;preparation_sha256=$ExpectedPreparationSha256;wrapper_sha256=$ExpectedWrapperSha256;launcher_sha256=$ExpectedLauncherSha256;authorization_sha256=(Get-FileHash -LiteralPath $authFull -Algorithm SHA256).Hash.ToLowerInvariant();production_admitted=$false}|ConvertTo-Json|Set-Content -LiteralPath ($prefix+'-preflight.json') -Encoding utf8
$arguments=@('-d','CFD-Foundation14','--exec','bash','/mnt/e/CFD-Learning-Lab/local/prospective-tracer-v001/run-openfoam-tracer-qualification-v004.sh',$Phase,$authLinux)
$child=Start-Process -FilePath 'wsl.exe' -ArgumentList $arguments -WorkingDirectory $private -WindowStyle Hidden -RedirectStandardOutput ($prefix+'.stdout.txt') -RedirectStandardError ($prefix+'.stderr.txt') -PassThru
$receipt=@{pid=$child.Id;phase=$Phase;time=(Get-Date -Format o);hidden=$true;processes=2;actual_steps_authorized=5;stdout=$prefix+'.stdout.txt';stderr=$prefix+'.stderr.txt';preflight=$prefix+'-preflight.json';production_admitted=$false;dedicated_distro_termination_deferred_until_native_commands_close_and_sync=$true}
$receipt|ConvertTo-Json|Set-Content -LiteralPath ($prefix+'-launch.json') -Encoding utf8
$receipt|ConvertTo-Json -Compress