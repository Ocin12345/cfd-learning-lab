# Dot-source this file in the shell used for this project.
$taskRoot = 'E:\CFD-Learning-Lab'
$env:TEMP = Join-Path $taskRoot 'temp'
$env:TMP = $env:TEMP
$env:TMPDIR = $env:TEMP
$env:UV_CACHE_DIR = Join-Path $env:TEMP 'uv-cache'
$env:UV_PYTHON_INSTALL_DIR = Join-Path $taskRoot 'environments\managed-python'
$env:PIP_CACHE_DIR = Join-Path $env:TEMP 'pip-cache'
$env:PYTHONPYCACHEPREFIX = Join-Path $env:TEMP 'pycache'
$env:PYTHONDONTWRITEBYTECODE = '1'
$env:PYTHONUTF8 = '1'
$env:MPLCONFIGDIR = Join-Path $env:TEMP 'matplotlib'
$env:PYFLUENT_LOGGING = 'OFF'
$env:FLUENT_EXECUTABLE = 'D:\ANSYS Inc\ANSYS Student\v261\fluent\ntbin\win64\fluent.exe'
$env:CFD_PROJECT_ROOT = $taskRoot
$env:CFD_OUTPUT_ROOT = Join-Path $taskRoot ('native\runs\' + (Get-Date -Format 'yyyy-MM-dd-HHmmss-fff'))
$env:PATH = (Join-Path $taskRoot 'environments\fluent\Scripts') + ';' + $env:PATH
foreach ($folder in @($env:TEMP,$env:UV_CACHE_DIR,$env:PIP_CACHE_DIR,$env:PYTHONPYCACHEPREFIX,$env:MPLCONFIGDIR)) {
    New-Item -ItemType Directory -Path $folder -Force | Out-Null
}
# Each dot-source selects a fresh output root. Keep the same shell while generating
# a mesh and solving it. No solver is started by this file.
