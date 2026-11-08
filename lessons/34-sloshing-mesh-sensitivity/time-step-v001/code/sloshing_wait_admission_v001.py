from pathlib import Path
import shutil,time
from sloshing_prelaunch_v001 import prelaunch,verify
from wing_fluent_v001 import available_memory
R=Path('E:/CFD-Learning-Lab')
def wait_admission(manifest,on_wait=lambda **v:None):
 while True:
  verify(manifest)
  if not Path('D:/ANSYS Inc/ANSYS Student/v261/fluent/ntbin/win64/fluent.exe').is_file():raise RuntimeError('Fluent executable missing')
  ram=available_memory();disk=shutil.disk_usage(R).free
  if ram<4*2**30 or disk<10*2**30:on_wait(reason='resources',RAM_bytes=ram,E_free_bytes=disk);time.sleep(15);continue
  try:return prelaunch(manifest)
  except AssertionError as e:
   if str(e) not in ['', 'Fluent/cortex already active or process gate failed','OpenFOAM/MPI active or native process check failed']:raise
   on_wait(reason='resource_race_or_process_gate',detail=str(e));time.sleep(15)
