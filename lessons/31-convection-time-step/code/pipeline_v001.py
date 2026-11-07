from pathlib import Path
import os,json,time,subprocess,shutil,traceback
from wing_fluent_v001 import available_memory
R=Path('E:/CFD-Learning-Lab');L=R/'local/natural-convection-time-step-v001';F=R/'native/runs/natural-convection-dt001-v001/fluent-v002';O=F.parent/'openfoam'
def mark(stage,**kw):(L/'pipeline-status.json').write_text(json.dumps(dict(stage=stage,**kw),indent=2)+'\n')
def run(name,args):
 mark(name)
 with (L/(name+'.log')).open('w',encoding='utf-8') as f:
  p=subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW,env=os.environ,stdin=subprocess.DEVNULL)
 assert p.returncode==0,(name,p.returncode)
try:
 mark('waiting_for_fluent_close')
 while True:
  status=json.loads((F/'status.json').read_text())
  assert status['stage']!='failed',status
  if status['stage']=='closed':break
  time.sleep(10)
 mark('resource_check')
 assert available_memory()>4*2**30 and shutil.disk_usage(R).free>10*2**30
 # The Fluent driver writes closed only after session.exit(). The shell also rejects an active foamRun.
 run('openfoam_run',['wsl.exe','-d','CFD-Foundation14','--','bash','/mnt/e/CFD-Learning-Lab/local/natural-convection-time-step-v001/run_openfoam_v001.sh'])
 run('openfoam_step_audit',[os.sys.executable,'-B',str(R/'local/natural-convection-production-v001/audit_openfoam_steps_v001.py'),str(O),'16','2000',str(L/'openfoam-steps.json')])
 run('openfoam_heat',['wsl.exe','-d','CFD-Foundation14','--','bash','/mnt/e/CFD-Learning-Lab/local/natural-convection-time-step-v001/postprocess_heat_v001.sh'])
 run('render',[os.sys.executable,'-B',str(L/'render_paired_v001.py')])
 run('comparison_audit',[os.sys.executable,'-B',str(L/'audit_comparison_v001.py')])
 mark('complete_ready_for_visual_review')
except Exception:
 mark('failed',error=traceback.format_exc());raise
