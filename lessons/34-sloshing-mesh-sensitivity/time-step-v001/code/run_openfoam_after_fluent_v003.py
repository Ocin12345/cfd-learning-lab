from pathlib import Path
from datetime import datetime
import os,time,json,hashlib,shutil,subprocess,traceback
from wing_fluent_v001 import available_memory
from sloshing_prelaunch_v001 import prelaunch,verify
from audit_sloshing_dt_qualification_v003 import audit
S=Path(__file__).parent;R=Path('E:/CFD-Learning-Lab');N=R/'native/runs/showcase-sloshing-time-step-v001';F=N/'fluent-qualification-v001';O=N/'openfoam-qualification-v001'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mark(stage,**v):
 p=S/'openfoam-qualification-v003-status.json';q=p.with_suffix('.tmp');q.write_text(json.dumps(dict(stage=stage,time=datetime.now().astimezone().isoformat(),pid=os.getpid(),**v),indent=2)+'\n');q.replace(p)
def main():
 assert not (S/'openfoam-qualification-v003-status.json').exists() and not (O/'solver.log').exists()
 f=json.loads((F/'qualification-audit-v001.json').read_text());assert f['numerical_pass'] and f['steps']==80
 while True:
  verify('qualification-inputs-v003.json')
  if available_memory()<4*2**30 or shutil.disk_usage(R).free<10*2**30:mark('waiting_resources');time.sleep(30);continue
  try:admission=prelaunch('qualification-inputs-v003.json')
  except AssertionError as e:
   if str(e) not in ['Fluent/cortex already active or process gate failed','OpenFOAM/MPI active or native process check failed']:raise
   mark('waiting_process_gate',reason=str(e));time.sleep(30);continue
  break
 for p,h in f['source_sha256'].items():assert sha(Path(p))==h,p
 verify('qualification-inputs-v003.json');mark('starting_first_OpenFOAM_80step_solve',admission=admission,accepted_Fluent_audit_sha256=sha(F/'qualification-audit-v001.json'))
 with (S/'openfoam-qualification-native-v003.log').open('xb') as log:p=subprocess.run(['wsl.exe','-d','CFD-Foundation14','-e','bash','/mnt/e/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/sloshing-time-step-v001/run_sloshing_dt_openfoam_qualification_v003.sh'],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
 assert p.returncode==0,p.returncode;verify('qualification-inputs-v003.json');mark('auditing_closed_OpenFOAM_qualification')
 a=audit('openfoam',O,O/'qualification-audit-v001.json');assert a['numerical_pass'];mark('paired_80step_qualifications_pass',states_each=81,end_time_s=.02,full_interval_not_launched=True,actual_visual_inspection_pending=True)
if __name__=='__main__':
 try:main()
 except Exception:mark('failed',error=traceback.format_exc());raise
