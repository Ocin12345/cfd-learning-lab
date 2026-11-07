from pathlib import Path
import os,json,time,subprocess,shutil,traceback
from wing_fluent_v001 import available_memory
R=Path('E:/CFD-Learning-Lab');L=R/'local/hot-patch-plume-production-v001';F=R/'native/runs/hot-patch-plume-20s-v001/fluent';O=F.parent/'openfoam'
def mark(stage,**kw):(L/'pipeline-status.json').write_text(json.dumps(dict(stage=stage,**kw),indent=2)+'\n')
def run(stage,args):
 mark(stage)
 with (L/(stage+'.log')).open('w',encoding='utf-8') as f:p=subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW,stdin=subprocess.DEVNULL,env=os.environ)
 assert p.returncode==0,(stage,p.returncode)
def sh(name):return ['wsl.exe','-d','CFD-Foundation14','--','bash','/mnt/e/CFD-Learning-Lab/local/hot-patch-plume-production-v001/'+name]
try:
 mark('waiting_for_fluent_close')
 while True:
  try:s=json.loads((F/'status.json').read_text())
  except json.JSONDecodeError:time.sleep(2);continue
  assert s['stage'] not in ['failed','stopped_after_save'],s
  if s['stage']=='closed':break
  time.sleep(10)
 fs=json.loads((F/'time-step-checks.json').read_text());assert len(fs)==2000 and all(x['passed'] for x in fs)
 assert available_memory()>4*2**30 and shutil.disk_usage(R).free>10*2**30
 run('openfoam_run',sh('run_openfoam_v001.sh'))
 run('openfoam_steps',[os.sys.executable,'-B',str(R/'local/natural-convection-production-v001/audit_openfoam_steps_v001.py'),str(O),'16','2000',str(L/'openfoam-steps.json')])
 run('openfoam_heat',sh('postprocess_openfoam_heat_v001.sh'))
 mark('both_runs_saved_ready_for_native_field_audit')
except Exception:
 mark('failed',error=traceback.format_exc());raise
