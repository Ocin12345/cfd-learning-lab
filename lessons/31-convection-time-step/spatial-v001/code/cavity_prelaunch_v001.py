from pathlib import Path
import hashlib,json,shutil,subprocess
from wing_fluent_v001 import available_memory
S=Path(__file__).resolve().parent;R=Path('E:/CFD-Learning-Lab')
def verify(name):
 a=json.loads((S/name).read_text())
 for p,h in a.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
 return len(a)
def prelaunch(name):
 count=verify(name)
 x=subprocess.run(['powershell','-NoProfile','-Command',"Get-CimInstance Win32_Process | Where-Object {$_.Name -in @('fluent.exe','cortex.exe')} | Select-Object ProcessId,CommandLine | ConvertTo-Json -Compress"],capture_output=True,text=True,creationflags=subprocess.CREATE_NO_WINDOW)
 assert x.returncode==0 and not x.stdout.strip(),'Fluent/cortex already active or process gate failed'
 x=subprocess.run(['wsl.exe','-d','CFD-Foundation14','-e','pgrep','-af','[f]oamRun|[s]impleFoam|[p]impleFoam|[m]pirun'],capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW);assert x.returncode==1,'OpenFOAM/MPI active or native process check failed'
 assert Path('D:/ANSYS Inc/ANSYS Student/v261/fluent/ntbin/win64/fluent.exe').is_file()
 verify(name);ram=available_memory();disk=shutil.disk_usage(R).free;assert ram>=4*2**30 and disk>=10*2**30
 return dict(verified_inputs=count,available_RAM_bytes=ram,free_E_bytes=disk,solver_ranks=2,one_solver=True)
