import subprocess,shutil,time,json,sys
sys.path.insert(0,"E:/CFD-Learning-Lab/repo/cfd-learning-lab/scripts")
from jet_energy_common_v001 import *
from wing_fluent_v001 import available_memory
MANIFEST=S/'full-interval-inputs-v001.json'
def verify():
 a=json.loads(MANIFEST.read_text())
 for p,h in a.items():
  if sha(p)!=h:raise RuntimeError('Frozen source changed: '+p)
 return len(a)
def prelaunch():
 count=verify()
 x=subprocess.run(['powershell','-NoProfile','-Command',"Get-CimInstance Win32_Process | Where-Object {$_.Name -in @('fluent.exe','cortex.exe')} | Select-Object ProcessId,CommandLine | ConvertTo-Json -Compress"],capture_output=True,text=True,creationflags=subprocess.CREATE_NO_WINDOW)
 if x.returncode or x.stdout.strip():raise BlockingIOError('Fluent/cortex active or process gate failed')
 x=subprocess.run(['wsl.exe','-d','CFD-Foundation14','-e','pgrep','-af','[f]oamRun|[s]impleFoam|[p]impleFoam|[m]pirun'],capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW)
 if x.returncode!=1 or x.stdout.strip():raise BlockingIOError('OpenFOAM/MPI active or process gate failed')
 if not Path('D:/ANSYS Inc/ANSYS Student/v261/fluent/ntbin/win64/fluent.exe').is_file():raise RuntimeError('Fluent executable absent')
 count=verify();ram=available_memory();disk=shutil.disk_usage(R).free
 if ram<4*2**30 or disk<10*2**30:raise BlockingIOError('RAM or E disk below threshold')
 return dict(verified_inputs=count,available_RAM_bytes=ram,E_free_bytes=disk,ranks=2,one_solver=True)
def wait_admission(notify=lambda **v:None):
 while True:
  try:return prelaunch()
  except BlockingIOError as e:notify(reason=str(e));time.sleep(15)
if __name__=='__main__':print(json.dumps(wait_admission()))
