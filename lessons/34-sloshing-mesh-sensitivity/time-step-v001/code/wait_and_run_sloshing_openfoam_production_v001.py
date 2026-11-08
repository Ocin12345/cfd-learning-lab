from pathlib import Path
from datetime import datetime
import os,json,hashlib,subprocess,shutil,time,traceback
from sloshing_dt_common_v001 import *
from sloshing_prelaunch_v001 import prelaunch,verify
from sloshing_full_audits_v002 import auditO
from wing_fluent_v001 import available_memory
F=N/'fluent-production-v001';O=N/'openfoam-production-v001'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mark(stage,**v):
 p=S/'openfoam-production-v001-status.json';q=p.with_suffix('.tmp');q.write_text(json.dumps(dict(stage=stage,time=datetime.now().astimezone().isoformat(),pid=os.getpid(),**v),indent=2)+'\n');q.replace(p)
def main():
 assert not (S/'openfoam-production-v001-status.json').exists() and not (O/'solver.log').exists()
 while True:
  p=S/'fluent-production-v001-status.json';a=json.loads(p.read_text()) if p.exists() else {}
  if a.get('stage')=='fluent_full_interval_pass':break
  if a.get('stage')=='failed':raise RuntimeError('Fluent full continuation failed; OpenFOAM full withheld')
  mark('waiting_for_accepted_full_Fluent',Fluent_stage=a.get('stage'),Fluent_step=a.get('step'));time.sleep(30)
 while True:
  verify('full-interval-inputs-v001.json')
  if available_memory()<4*2**30 or shutil.disk_usage(R).free<10*2**30:mark('waiting_resources_after_Fluent');time.sleep(30);continue
  try:admission=prelaunch('full-interval-inputs-v001.json')
  except AssertionError as e:
   if str(e) not in ['Fluent/cortex already active or process gate failed','OpenFOAM/MPI active or native process check failed']:raise
   mark('waiting_process_cleanup_gate',reason=str(e));time.sleep(30);continue
  break
 a=json.loads((F/'production-audit-v001.json').read_text());assert a['numerical_pass'] and a['merged_steps']==4800
 for p,h in a['source_sha256'].items():assert sha(Path(p))==h,p
 verify('full-interval-inputs-v001.json');assert available_memory()>=4*2**30 and shutil.disk_usage(R).free>=10*2**30;mark('starting_exact_OpenFOAM_continuation',admission=admission,accepted_Fluent_audit_sha256=sha(F/'production-audit-v001.json'))
 with (S/'openfoam-production-native-v001.log').open('xb') as log:p=subprocess.run(['wsl.exe','-d','CFD-Foundation14','-e','bash','/mnt/e/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/sloshing-time-step-v001/run_sloshing_openfoam_production_v001.sh'],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
 assert p.returncode==0,p.returncode;verify('full-interval-inputs-v001.json');mark('auditing_closed_OpenFOAM_full_interval');a=auditO(O);assert a['numerical_pass'];mark('openfoam_full_interval_pass',new_steps=4720,combined_steps=4800,paired_media_not_yet_exported=True)
if __name__=='__main__':
 try:main()
 except Exception:mark('failed',error=traceback.format_exc());raise
