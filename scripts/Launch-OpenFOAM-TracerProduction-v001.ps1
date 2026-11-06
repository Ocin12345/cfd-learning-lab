# SOURCE ONLY until exact root admission. Hidden795-step accepted-pulse continuation.
param(
 [Parameter(Mandatory=$true)][string]$Authorization,
 [Parameter(Mandatory=$true)][string]$ExpectedPreparationSha256,
 [Parameter(Mandatory=$true)][string]$ExpectedWrapperSha256,
 [Parameter(Mandatory=$true)][string]$ExpectedLauncherSha256
)
$ErrorActionPreference='Stop'
. 'E:\CFD-Learning-Lab\local\Use-Project.ps1'
$private='E:\CFD-Learning-Lab\local\prospective-tracer-v001'
$prep=Join-Path $private 'prepare_openfoam_tracer_production_v001.py'
$wrapper=Join-Path $private 'run-openfoam-tracer-production-v001.sh'
$python='E:\CFD-Learning-Lab\environments\fluent\Scripts\python.exe'
$authFull=(Resolve-Path -LiteralPath $Authorization).Path
if(!$authFull.StartsWith($private+'\',[StringComparison]::OrdinalIgnoreCase)){throw 'Authorization must remain in exact candidate private folder'}
$auth=Get-Content -LiteralPath $authFull -Raw | ConvertFrom-Json
if($auth.approved -ne $true -or $auth.approved_by -cne 'root' -or $auth.action -cne 'run_tracer_openfoam_production_v001' -or $auth.additional_steps -ne 795 -or $auth.pilot_steps -ne 5 -or $auth.total_steps -ne 800 -or $auth.total_states -ne 801 -or $auth.start_time_s -ne 0.125 -or $auth.final_time_s -ne 20 -or $auth.deltaT_s -ne 0.025 -or $auth.processes -ne 2 -or $auth.fluent_sessions_closed -ne $true -or $auth.fluent_qualification_accepted_by_root -ne $true -or $auth.openfoam_control_accepted_by_root -ne $true -or $auth.openfoam_pulse_accepted_by_root -ne $true){throw 'Exact root admission after both codes pilots and separate new OF control/pulse required'}
foreach($pair in @(@($prep,$ExpectedPreparationSha256),@($wrapper,$ExpectedWrapperSha256),@($PSCommandPath,$ExpectedLauncherSha256))){
 if($pair[1] -cnotmatch '^[a-f0-9]{64}$' -or (Get-FileHash -LiteralPath $pair[0] -Algorithm SHA256).Hash.ToLowerInvariant() -cne $pair[1]){throw ('Reviewed source changed: '+$pair[0])}
}
if($auth.preparation_sha256 -cne $ExpectedPreparationSha256 -or $auth.wrapper_sha256 -cne $ExpectedWrapperSha256 -or $auth.launcher_sha256 -cne $ExpectedLauncherSha256){throw 'Authorization must bind all3 exact reviewed sources'}
if(@(Get-Process fluent*,cortex*,fl_mpi*,fl2d*,fl3d* -ErrorAction SilentlyContinue).Count){throw 'Fluent must be closed before any WSL/native preflight'}
$ram=[int64](Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory*1024
$disk=[int64](Get-PSDrive E).Free
if($ram -lt 4GB -or $disk -lt 10GB){throw ('Resource admission held: '+[math]::Round($ram/1GB,3)+' GiB RAM; '+[math]::Round($disk/1GB,3)+' GiB E:')}
$casePath='E:\CFD-Learning-Lab\native\runs\2026-11-02-tracer-tjunction-v001\openfoam-production-v001'
if(Test-Path -LiteralPath $casePath){throw 'Fresh production case required; preserve existing attempts'}
& $python -B $prep --authorization $authFull --verify-sources
if($LASTEXITCODE -ne 0){throw 'Accepted pilot/source/auth preflight failed before WSL start'}
# Source audits may take time: read the actual RAM/disk and Fluent processes again immediately before launch.
$ram=[int64](Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory*1024
$disk=[int64](Get-PSDrive E).Free
if($ram -lt 4GB -or $disk -lt 10GB -or @(Get-Process fluent*,cortex*,fl_mpi*,fl2d*,fl3d* -ErrorAction SilentlyContinue).Count){throw 'Fresh final prelaunch resource/sequential-solver admission held'}
$stamp=Get-Date -Format 'yyyyMMdd-HHmmss-fff'
$prefix=Join-Path $private ('openfoam-production-v001-'+$stamp)
$authLinux='/mnt/e/'+$authFull.Substring(3).Replace('\','/')
@{time=(Get-Date -Format o);available_ram_bytes=$ram;ram_floor_bytes=4GB;free_E_bytes=$disk;disk_floor_bytes=10GB;processes=2;fluent_active=$false;pilot_steps=5;additional_steps=795;total_steps=800;start_time_s=0.125;final_time_s=20;production_targets_sha256=$auth.production_targets_sha256;preparation_sha256=$ExpectedPreparationSha256;wrapper_sha256=$ExpectedWrapperSha256;launcher_sha256=$ExpectedLauncherSha256;authorization_sha256=(Get-FileHash -LiteralPath $authFull -Algorithm SHA256).Hash.ToLowerInvariant();full_signed_scalar_budget_analysis_pending=$true}|ConvertTo-Json|Set-Content -LiteralPath ($prefix+'-preflight.json') -Encoding utf8
$arguments=@('-d','CFD-Foundation14','--exec','bash','/mnt/e/CFD-Learning-Lab/local/prospective-tracer-v001/run-openfoam-tracer-production-v001.sh',$authLinux)
$child=Start-Process -FilePath 'wsl.exe' -ArgumentList $arguments -WorkingDirectory $private -WindowStyle Hidden -RedirectStandardOutput ($prefix+'.stdout.txt') -RedirectStandardError ($prefix+'.stderr.txt') -PassThru
$receipt=@{pid=$child.Id;time=(Get-Date -Format o);hidden=$true;processes=2;pilot_steps=5;additional_steps_authorized=795;total_steps=800;stdout=$prefix+'.stdout.txt';stderr=$prefix+'.stderr.txt';preflight=$prefix+'-preflight.json';full_signed_scalar_budget_analysis_pending=$true;dedicated_distro_termination_deferred_until_native_commands_close_and_sync=$true}
$receipt|ConvertTo-Json|Set-Content -LiteralPath ($prefix+'-launch.json') -Encoding utf8
$receipt|ConvertTo-Json -Compress
