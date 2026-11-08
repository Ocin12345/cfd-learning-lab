import sys,subprocess,time,shutil,json
sys.path.insert(0,'E:/CFD-Learning-Lab/repo/cfd-learning-lab/scripts')
from wing_repair_common_v001 import *
from wing_fluent_v001 import available_memory
MANIFEST=W/'retained-probe-inputs-v001.json'
def verify():
 d=json.loads(MANIFEST.read_text())
 for p,h in d.items():
  if sha(p)!=h:raise RuntimeError('Changed frozen wing source/parent: '+p)
 return len(d)
def prelaunch():
 n=verify();exe=Path('D:/ANSYS Inc/ANSYS Student/v261/fluent/ntbin/win64/fluent.exe')
 if not exe.is_file():raise RuntimeError('Fluent executable absent')
 p=subprocess.run(['powershell.exe','-NoProfile','-Command',"Get-CimInstance Win32_Process | Where-Object {$_.Name -in @('fluent.exe','cortex.exe')} | Select-Object ProcessId,CreationDate,Name | ConvertTo-Json -Compress"],capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW)
 if p.returncode or p.stdout.strip():raise BlockingIOError('Windows Fluent/cortex conflict or failed process gate')
 p=subprocess.run(['wsl.exe','-d','CFD-Foundation14','-e','pgrep','-af','[f]oamRun|[s]impleFoam|[p]impleFoam|[m]pirun|[r]econstructPar'],capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW)
 if p.returncode!=1 or p.stdout.strip():raise BlockingIOError('WSL foam/MPI/reconstruction conflict or failed process gate')
 n=verify();ram=available_memory();disk=shutil.disk_usage(R).free
 if ram<4*2**30 or disk<10*2**30:raise BlockingIOError('RAM/Efree below4GiB/10GiB')
 return dict(verified_inputs=n,RAM_bytes=ram,E_free_bytes=disk,ranks=2,one_solver=True)
def wait_admission(notify=lambda **v:None):
 while True:
  try:return prelaunch()
  except BlockingIOError as e:notify(reason=str(e));time.sleep(15)
